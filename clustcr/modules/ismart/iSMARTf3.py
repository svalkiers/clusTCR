#! usr/bin/python3
# -*- coding: utf-8 -*-
## Ultra-fast pairwise alignment algorithm to analyze up to 10^8 CDR3s

import sys, os, re, resource
import numpy as np
from substitution_matrices import BLOSUM62
import time
from operator import itemgetter
from itertools import chain
from random import shuffle
from optparse import OptionParser
from collections import Counter

AAstring='ACDEFGHIKLMNPQRSTVWY'
AAstringList=list(AAstring)
cur_dir=os.path.dirname(os.path.realpath(__file__))+'/'

class CDR3:
    def __init__(self, s, sID, KS, st, ed):
        ## initialize with an input sequence
        ## s: input CDR3 sequence starting from C and ending with the first F in FGXG
        ## sID: unique identifier (increasing integers) given to each CDR3 sequence. Even identical CDR3s should have distinct sIDs
        ## KS: Kmer size
        ## st: the first 0:(st-1) amino acids will not be included in K-merization
        ## ed: the last L-ed amino acids will be skipped
        self.s=s
        self.ID=sID
        L=len(s)
        self.L=L
        sub_s=s[st: (L-ed)]
        Ls=len(sub_s)
        Kmer=[sub_s[x:(x+KS)] for x in range(0,Ls-KS+1)]
        self.Kmer=Kmer

class KmerSet:
    ## Kmer set for fast read searching based on mismatch-allowed Kmer index
    def __init__(self, Seqs, sIDs, KS, st, ed):
        ## initialize with a list of CDR3s, parse each CDR3 into Kmers, build Kmer-sID dictionary
        ## Seqs and sIDs must have the same length
        if len(Seqs) != len(sIDs):
            raise "Sequence and ID lists have different length. Please check input."
        KmerDict={}
        N=len(Seqs)
        self.N=N
        CDR3Dict={}
        LLs=[]
        for ii in range(0,N):
            s=Seqs[ii]
            sID=sIDs[ii]
            cc=CDR3(s,sID,KS,st,ed)
            CDR3Dict[cc.ID]=cc.Kmer
            KK=cc.Kmer
            LLs.append(cc.L)
            for kk in KK:
                if kk not in KmerDict:
                    KmerDict[kk]=[sID]
                else:
                    KmerDict[kk].append(sID)
        self.KD=KmerDict
        self.KS=KS
        self.CD=CDR3Dict
        self.LL=LLs
    def FindKmerNeighbor(self,kk):
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            kk_pre=[kk[0:jj]]*20
            kk_suf=[kk[(jj+1):KS]]*20
            kkn=list(zip(kk_pre,AAstringList,kk_suf))
            KS_n1+=[''.join(list(x)) for x in kkn]
        return KS_n1
    def FindKmerNeighbor2(self,kk):
        ## KS≥6, allowing 2 mismatches. CDR3 length must be ≥ 10
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            for ii in range(KS):
                if ii<=jj:
                    continue
                kk_pre=[kk[0:jj]]*20
                kk_mid=[kk[(jj+1):ii]]*20
                kk_suf=[kk[(ii+1):KS]]*400
                kkn=list(zip(kk_pre,AAstringList,kk_mid))
                kkn=[''.join(list(x)) for x in kkn]
                kkn=[[x]*20 for x in kkn]
                kkn=list(chain(*kkn))
                kkn2=list(zip(kkn, AAstringList*20, kk_suf))
                kkn2=[''.join(list(x)) for x in kkn2]
                KS_n1+=kkn2
        return KS_n1
    def KmerIndex(self):
        ## For each K-mer, find its nearest neighbor with 1 character mismatch
        KKs=list(self.KD.keys())
        KS=self.KS
        KKs_set=set(KKs)
        Skk='_'.join(KKs)
        KI_Dict={}
        for kk in KKs:
##            kk_neighbor=[]
##            for jj in range(KS):
##                kk_pre=kk[0:jj]
##                kk_suf=kk[(jj+1):KS]
##                pat=kk_pre+'['+AAstring+']{1}'+kk_suf
##                p=re.compile(pat)
##                mm=[m.group() for m in p.finditer(Skk)]
##                kk_neighbor+=mm
            KS_n=set(self.FindKmerNeighbor(kk))
            kk_neighbor = KS_n & KKs_set
            KI_Dict[kk]=list(kk_neighbor)
        return KI_Dict
    def updateKD(self, KI):
        ## group sequences sharing motifs with 1-2 mismatches
        KD=self.KD
        KDnew={}
        for kk in KD:
            kkm=KI[kk]
            vvL=itemgetter(*kkm)(KD)
            if isinstance(vvL[0],list):
                vvL=list(chain(*vvL))
            KDnew[kk]=vvL
        return KDnew

def GenerateMotifGraph(mD,seqs,seqID):
    SeqShareGraph={}
    mDL={}
    for kk in mD:
        vv=mD[kk]
        LL=[]
        for v in vv:
            LL.append(len(seqs[v]))
        mDL[kk]=LL
    for kk in mD:
        vv=mD[kk]
        LL=mDL[kk]
        nv=len(vv)
        for ii in range(0,nv):
            id_1=vv[ii]
            L1=LL[ii]
            for jj in range(ii,nv):
                if jj==ii:
                    continue
                id_2=vv[jj]
                L2=LL[jj]
                if L2 != L1:
                    continue
                if id_1 not in SeqShareGraph:
                    SeqShareGraph[id_1]=[id_2]
                elif id_2 not in SeqShareGraph[id_1]:
                    SeqShareGraph[id_1].append(id_2)
                if id_2 not in SeqShareGraph:
                    SeqShareGraph[id_2]=[id_1]
                elif id_1 not in SeqShareGraph[id_2]:
                    SeqShareGraph[id_2].append(id_1)
    return SeqShareGraph

def generateSSG(Kset, CDR3s, k_thr=2):
    KD=Kset.KD
    KI=Kset.KmerIndex()
    KDnew=Kset.updateKD(KI)
    CD=Kset.CD
    LL=np.array(Kset.LL)
    SSG={}
    for kk in CD:
        vv=itemgetter(*CD[kk])(KDnew)
        if isinstance(vv[0],list):
            vv=list(chain(*vv))
        vv1=[]
        c=Counter(vv)
        for k in c:
            if c[k]>=k_thr:
                vv1.append(k)
        vv1=np.array(vv1)
        if len(vv1)==0:
            continue
        cdr3=CDR3s[kk]
        L0=len(cdr3)
        idx=np.where(LL[vv1]==L0)[0]
        if len(idx)==0:
            continue
        vvs=list(vv1[idx])
        vvs.remove(kk)
        if len(vvs)>0:
            SSG[kk]=vvs
    return SSG

def UpdateSSG(SSG, seqs, Vgenes, Vscore={}, UseV=True, gap=-6, gapn=1, cutoff=7.5):
    SSGnew={}
    count=0
    t1=time.time()
    N=len(list(chain(*list(SSG.values()))))
    print("Number of pairs to be processed: %d" %N)
    for kk in SSG:
        s1=seqs[kk]
        V1=Vgenes[kk]
        VV=SSG[kk]
        for vv in VV:
            s2=seqs[vv]
            V2=Vgenes[vv]
            score=falign(s1, s2, V1, V2, st=2, VScore=Vscore, UseV=UseV, gap=-6, gapn=1)
            count+=1
            if count % 1000000 ==0:
                t2=time.time()
                print("Processed %d pairs. Elapsed time %f" %(count, t2-t1))
            if score>=cutoff:
                if kk not in SSGnew:
                    SSGnew[kk]=[vv]
                else:
                    SSGnew[kk].append(vv)
    return SSGnew
            
##def IdentifyMotifCluster0(SSG):
##    ## Input SeqShareGraph dictionary representation of sparse matrix
##    POS=SSG.keys()
##    NP=len(POS)
##    ClusterList=[]
##    tmpL=list(chain(*ClusterList))
##    count=0
##    for ii in POS:
##        if ii not in tmpL:
###            STACK=LoadComm([],ii)
##            STACK=dfs(SSG,ii)
##            ClusterList.append(list(STACK))
##            tmpL=list(chain(*ClusterList))
##            count+=1
##            if count % 200 ==0:
##                print ("    Solved %d clusters" %(count))
##    return ClusterList

def IdentifyMotifCluster(SSG):
    ## Input SeqShareGraph dictionary representation of sparse matrix
    POS=set(SSG.keys())
    NP=len(POS)
    ClusterList=[]
    tmpL=set(chain(*ClusterList))
    count=0
    while 1:
            xx=POS ^ tmpL
            if len(xx)==0:
                break
            for ii in xx:
#            STACK=LoadComm([],ii)
                STACK=dfs(SSG,ii)
                tmpL = tmpL | STACK
                ClusterList.append(list(STACK))
#                tmpL=set(chain(*ClusterList))
                count+=1
                if count % 200 ==0:
                    print ("    Solved %d clusters" %(count))
                break
    return ClusterList

def dfs(graph, start):
    '''
    Non-resursive depth first search
    '''
    visited = set()
    stack = [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(set(graph[vertex]) - visited)
    
    return visited

def ParseInput(filename, KS, UseV=True):
    ## Input file format (in column order): CDR3, Variable Gene, Frequency, Other Info (tab deliminated)
    g=open(filename)
    DataDict={}
    CDR3s=[]
    Vgene=[]
    count=0
    flag=0
    Header=''
    while 1:
        line=g.readline()
        if len(line)==0:
            break
        ww=line.strip().split('\t')
        if ww[0].startswith('C') and ww[0].endswith('F'):
            s=ww[0]
            if len(s) < KS+4:
                continue
            DataDict[count]=ww
            count+=1
            CDR3s.append(ww[0])
            if UseV:
                Vgene.append(ww[1])
        elif flag==0:
            Header=line
            flag=1
    sIDs=list(DataDict.keys())
    return CDR3s, Vgene, sIDs, DataDict, Header

def InsertGap(Seq,n):
    ## Insert n gaps to Seq; n<=2
    if n==0:
        return [Seq]
    ns=len(Seq)
    SeqList=[]
    if(n==1):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            SeqList.append(SeqNew)
    if(n==2):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            for jj in range(0,ns+2):
                SeqNew0=SeqNew[0:jj]+'-'+SeqNew[jj:]
                SeqList.append(SeqNew0)
    return SeqList

def ParseFa(fname):
    InputStr=open(fname).readlines()
    FaDict={}
    seq=''
    for line in InputStr:
        if line.startswith('>'):
            if len(seq)>0:
                FaDict[seqHead]=seq
                seq=''
            seqHead=line.strip()
        else:
            seq+=line.strip()
    if seqHead not in FaDict:
        FaDict[seqHead]=seq
    return FaDict

def PreCalculateVgeneDist(VgeneFa="Imgt_Human_TRBV.fasta"):
    ## Only run one time if needed
    FaDict=ParseFa(cur_dir+VgeneFa)
    VScore={}
    CDR1Dict={}
    CDR2Dict={}
    for kk in FaDict:
        if '|' in kk:
            VV=kk.split('|')[1]
        else:
            VV=kk[1:]
        CDR1Dict[VV]=FaDict[kk][26:37]  ## Imgt CDR1: 27 - 38
        CDR2Dict[VV]=FaDict[kk][55:64]  ## Imgt CDR2: 56 - 65
    Vkeys=list(CDR1Dict.keys())
    nn=len(Vkeys)
    for ii in range(0,nn):
        V1=Vkeys[ii]
        s1_CDR1=CDR1Dict[V1]
        s1_CDR2=CDR2Dict[V1]
        for jj in range(ii,nn):
            V2=Vkeys[jj]
            s2_CDR1=CDR1Dict[V2]
            s2_CDR2=CDR2Dict[V2]
            score1=SeqComparison(s1_CDR1,s2_CDR1)
            score2=SeqComparison(s2_CDR2,s2_CDR2)
            #print score1+score2
            VScore[(V1,V2)]=score1+score2
    gg=open('VgeneScores.txt','w')
    for kk in VScore:
        vv=VScore[kk]
        line=kk[0]+'\t'+kk[1]+'\t'+str(vv)+'\n'
        gg.write(line)
    gg.close()

def SeqComparison(s1,s2,gap=-6):
    n=len(s1)
    CorList=[]
    score=0
    for kk in range(0,n):
        aa=s1[kk]
        bb=s2[kk]
        if aa in ['.','-','*'] or bb in ['.','-','*']:
            if aa!=bb:
                score += gap
            continue
        if aa==bb:
            score += min(4,BLOSUM62[(aa,aa)])
            continue
        KEY=(aa,bb)
        if KEY not in BLOSUM62:
            KEY=(bb,aa)
        if KEY not in BLOSUM62:
            raise "Non-standard amino acid coding!"
        score+=BLOSUM62[KEY]
    return score

def NHLocalAlignment(Seq1,Seq2,gap_thr=1,gap=-6):
    n1=len(Seq1)
    n2=len(Seq2)
    if n1<n2:
        Seq=Seq1
        Seq1=Seq2
        Seq2=Seq
        nn=n2-n1
    else:
        nn=n1-n2
    if nn>gap_thr:
        return -1
    SeqList1=[Seq1]
    SeqList2=InsertGap(Seq2,nn)
    alns=[]
    SCOREList=[]
    for s1 in SeqList1:
        for s2 in SeqList2:
                SCOREList.append(SeqComparison(s1,s2,gap))
    maxS=max(SCOREList)
    return maxS

def fun_map(p,f):
    ## Fake function for passing multiple arguments to Pool.map()
    return f(*p)

def falign(s1, s2, V1, V2 ,st,VScore={}, UseV=True, gapn=1, gap=-6):
    mid1=s1[st:-2]
    mid2=s2[st:-2]
    if UseV:
        if V2==V1:
            V_score=4
        else:
            Vkey=(V1,V2)
            if Vkey not in VScore:
                Vkey=(V2,V1)
            if Vkey not in VScore:
                #print("V gene not found!")
                return 0
            else:
                V_score=VScore[Vkey]/20.0
    else:
        V_score=4.0
    aln=NHLocalAlignment(mid1,mid2,gapn,gap)
    score=aln/float(max(len(mid1),len(mid2)))+V_score
    return score

def PWalign(Seqs,ID,Vgene={}, VScore={}, gap=-6,gapn=1,UseV=True,cutoff=7,Nthread=1):
    ## Wrapper function
    ns=len(Seqs)
    if ns != len(Vgene):
        if len(Vgene)==0:
            Vgene=['']*ns
        else:
            raise "Incompatible variable gene number!"
    z=sorted(zip(Seqs,Vgene,ID),key=lambda pair:len(pair[0]))
    Seqs=[x for x,y,t in z]
    Vgene=[x for y,x,t in z]
    ID=[x for t,y,x in z]
    del z
    PWscore={}
    st=4
    if not UseV:
        st=2
    if Nthread==1:
        for ii in range(0,ns):
            V1=Vgene[ii]
            for jj in range(ii,ns):
                if ii==jj:
                    continue
                V2=Vgene[jj]
                mid1=Seqs[ii][st:-2]
                mid2=Seqs[jj][st:-2]
                if UseV:
                    if V2==V1:
                        V_score=4
                    else:
                        Vkey=(V1,V2)
                        if Vkey not in VScore:
                            Vkey=(V2,V1)
                        if Vkey not in VScore:
                            #print("V gene not found!")
                            continue
                        else:
                            V_score=VScore[Vkey]/20.0  ## Take the floor of the float number
                else:
                    V_score=4.0
                aln=NHLocalAlignment(mid1,mid2,gapn,gap)
                #print aln
    #            J_score=NHLocalAlignment(Jend1,Jend2,gap=False)[0]
                score=aln/float(max(len(mid1),len(mid2)))+V_score
                if score>=cutoff:
                    PWscore[(ii,jj)]=1
    else:
        # Multi-thread processing
        p=Pool(Nthread)
        XX=[]
        for ii in range(0,ns):
            for jj in range(ii,ns):
                if ii==jj:
                    continue
                else:
                    XX.append([ii,jj])
        para= []
        for xx in XX:
            para.append((xx,st,VScore,Seqs, Vgene, UseV, gapn, gap))
        pl_out=p.map(partial(fun_map,f=falign),para)
        p.close()
        p.join()
        ## End multiple processing
        for kk in range(0,len(XX)):
            score=pl_out[kk]
            if score>=cutoff:
                PWscore[(XX[kk][0],XX[kk][1])]=1
    return (PWscore,Seqs,Vgene,ID)  

def IdentifyCDR3Clusters(PWscore,cutoff=7):
    POS=np.array(list(PWscore.keys()))[np.where(np.array(list(PWscore.values()))==1)]
    if len(POS)<=0:
        #print "Too few clustered CDR3s! Please check your repertoire data."
        return []
    POS=list(POS)
    POS=np.array([list(map(lambda x:x[0],POS)),list(map(lambda x:x[1],POS))])
    uniquePos=list(set(list(POS[0])+list(POS[1])))
    ClusterList=[]
    tmpL=list(chain(*ClusterList))
    def LoadComm(STACK,cur_ii):
        if cur_ii in STACK:
                return
        else:
                STACK.append(cur_ii)
                vv=list(POS[1][np.where(POS[0]==cur_ii)])+list(POS[0][np.where(POS[1]==cur_ii)])
                for v in vv:
                        LoadComm(STACK,v)
        return STACK
    for ii in uniquePos:
        if ii in tmpL:
            continue
        else:
            STACK=LoadComm([],ii)
            ClusterList.append(STACK)
            tmpL=list(chain(*ClusterList))
    return ClusterList

def runMotifClustering(f, KS, st, ed, THR, outDIR=cur_dir, Vscore={},UseV=True, k_thr=1, gap=-6, gapn=1, cutoff=7, Nthread=1):
    ## Wrapper function to perform k-mer based partition
    t1=time.time()
    CDR3s, Vgenes, sIDs, DataDict, Header = ParseInput(f, KS=KS, UseV=UseV)
    if not UseV:
        Vgenes=['TRBV2']*len(CDR3s)
    print("Collecting data and building K-mer index")
    Kset=KmerSet(CDR3s, sIDs, KS, st, ed)
    SSG=generateSSG(Kset, CDR3s, k_thr=k_thr)
    t2=time.time()
    print("Done! Time Elapsed %f" %(t2-t1))
    print("Performing pairwise alignment")
    SSGnew=UpdateSSG(SSG, CDR3s, Vgenes)
    print("Done!")
    print("Dividing CDR3s into clusters")
    CLall=IdentifyMotifCluster(SSGnew)
    t2=time.time()
    print("Done! Time Elapsed %f" %(t2-t1))
    f_name=f.split('/')
    f_name=f_name[len(f_name)-1]
    f_prefix=re.sub('\\.[txcsv]+','',f_name)
    f_out=outDIR+'/'+f_prefix+'_clustered_v3.txt'
    g=open(f_out,'w')
    InfoLine='#'+f+'|KS='+str(KS)+'|k_thr='+str(k_thr)+'|st='+str(st)+'|ed='+str(ed)+'|THR='+str(THR)+'|UseV='+str(int(UseV))+'|gap='+str(gap)+'|gapn='+str(gapn)+'|cutoff='+str(cutoff)
    g.write(InfoLine+'\n')
    if len(Header)>0:
       g.write(Header+'\t'+'Group'+'\n')
    for ii in range(len(CLall)):
        v=CLall[ii]
        for jj in v:
            info=DataDict[jj]
            line='\t'.join(info)+'\t'+str(ii)+'\n'
            g.write(line)
    g.close()    
    
def runClustering(f, KS, st, ed, THR, outDIR=cur_dir,Vscore={}, UseV=True, gap=-6, gapn=1, cutoff=7, Nthread=1):
    ## Wrapper function to perform k-mer based partition, pairwise alignment and clustering
    ## Input filename
    t1=time.time()
    CDR3s, Vgenes, sIDs, DataDict, Header = ParseInput(f, UseV=UseV)
    if not UseV:
        Vgenes=['TRBV2']*len(CDR3s)
    print("Creating sequence partition.")
    FP, LP = RecursivePartition(CDR3s, sIDs, KS=KS, st=st, ed=ed, THR=THR)
    t2=time.time()
    print("Done! Elapsed time %f" %(t2-t1))
    CLall=[]
    for fp in FP+LP:
        ID=fp
        if len(ID)>=THR:
            print("..processing large cluster %d" %len(ID))
        Seqs=list(np.array(CDR3s)[np.array(ID)])
        VG=list(np.array(Vgenes)[np.array(ID)])
        (PWs, Seqs, VG, ID)=PWalign(Seqs, ID, VG, Vscore, UseV= UseV, gap=gap, gapn=gapn, cutoff=cutoff, Nthread=Nthread)
        CL=IdentifyCDR3Clusters(PWs)
        for cl in CL:
            CLall.append(list(np.array(ID)[np.array(cl)]))
        if len(ID)>=THR/2:
            t2=time.time()
            print("...Elapsed time %f" %(t2-t1))
    f_name=f.split('/')
    f_name=f_name[len(f_name)-1]
    f_prefix=re.sub('\\.[txcsv]+','',f_name)
    f_out=outDIR+'/'+f_prefix+'_clustered_v3.txt'
    g=open(f_out,'w')
    InfoLine='#'+f+'|KS='+str(KS)+'|st='+str(st)+'|ed='+str(ed)+'|THR='+str(THR)+'|UseV='+str(int(UseV))+'|gap='+str(gap)+'|gapn='+str(gapn)+'|cutoff='+str(cutoff)
    g.write(InfoLine+'\n')
    if len(Header)>0:
       g.write(Header+'\t'+'Group'+'\n')
    for ii in range(len(CLall)):
        v=CLall[ii]
        for jj in v:
            info=DataDict[jj]
            line='\t'.join(info)+'\t'+str(ii)+'\n'
            g.write(line)
    g.close()

def CommandLineParser():
    parser=OptionParser()
    print ('''
iSMARTf is intended to perform pairwise CDR3 alignment for large volume (10^6-10^8) of sequneces. It implements
k-mer based clustering to recursively reduce the search space, and report sequences with high similarity.
Currently, gap is not supported. 

Input columns:
1. CDR3 amino acid sequence (Starting from C, ending with the first F/L in motif [FL]G.G)
2. Variable gene name in Imgt format: TRBVXX-XX*XX
3. Joining gene name (optional)
4. Frequency (optional)
5. Other information (optional)
''')
    parser.add_option("-d","--directory",dest="Directory",help="Input repertoire sequencing file directory. Please make sure that all the files in the directory are input files.",default="")
    parser.add_option("-f","--file",dest="File",default='',help="Input single file of CDR3 sequences for grouping")
    parser.add_option("-F","--fileList",dest="files",default='',help='Alternative input: a file containing the full path to all the files. If given, overwrite -d and -f option')
    parser.add_option("-T","--KmerThreshold",dest='THR',default=1000, help='Maximum k-mer cluster size.')
    parser.add_option("-k","--KmerNum", dest='kt', default=1, help="Number of k-mers to cross-index similar CDR3s.")
    parser.add_option("-K","--KmerSize",dest='KS',default=6,help='Kmer length')
    parser.add_option("-t","--threshold",dest="thr",default=7.5,help="Threshold for calling similar CDR3 groups. The higher the more specific.")
    parser.add_option("-o","--output",dest="OutDir",default='./',help="Output directory for intermediate and final outputs.")
    parser.add_option("-g","--GapPenalty",dest="Gap",default= -6,help="Gap penalty,default= -6")
    parser.add_option("-n","--GapNumber",dest="GapN",default=1,help="Maximum number of gaps allowed when performing alignment. Max=1, default=1")
    parser.add_option("-V","--VariableGeneFa",dest="VFa",default="Imgt_Human_TRBV.fasta",help="IMGT Human beta variable gene sequences")
    parser.add_option("-v","--VariableGene",dest="V",default=True,action="store_false",help="If False, iSMART will omit variable gene information and use CDR3 sequences only. This will yield reduced specificity. The cut-off will automatically become the current value-4.0")
    parser.add_option("-N","--NumberOfThreads",dest="NN",default=1,help="Number of threads for multiple processing. Not working so well.")
#    parser.add_option("-D","--UseDiAAmat",dest="Di",default=False,action="store_true",help="If True, iSMART will use a predefined di-amino acid substitution matrix in sequence comparison.")
    return parser.parse_args()

def main():
    (opt,_)=CommandLineParser()
    FileDir=opt.Directory
    if len(FileDir)>0:
            files=os.listdir(FileDir)
            files0=[]
            for ff in files:
                    ff=FileDir+'/'+ff
                    files0.append(ff)
            files=files0
    else:
            files=[]
    File=opt.File
    if len(File)>0:
            files=[File]
    FileList=opt.files
    if len(FileList)>0:
            files=[]
            fL=open(FileList)
            for ff in fL.readlines():
                    files.append(ff.strip())
    VFa=opt.VFa
    PreCalculateVgeneDist(VFa)
    vf=open('./VgeneScores.txt')  ## Use tcrDist's Vgene 80-score calculation
    VScore={}
    VV=opt.V
    if VV:
        while 1:
            line=vf.readline()
            if len(line)==0:
                break
            ww=line.strip().split('\t')
            VScore[(ww[0],ww[1])]=int(ww[2])
    Gap=int(opt.Gap)
    Gapn=int(opt.GapN)
    cutoff=float(opt.thr)
    OutDir=opt.OutDir
    st=3
    ed=1
    KS=int(opt.KS)
    THR=float(opt.THR)
    NT=int(opt.NN)
    kt=int(opt.kt)
    for ff in files:
        print("Processing %s" %ff)
        runMotifClustering(ff,KS=KS, st=st, ed=ed, THR=THR, outDIR=OutDir, Vscore=VScore,k_thr=kt, UseV=VV, gap=Gap, gapn=Gapn, cutoff=cutoff, Nthread=NT)

if __name__ == "__main__":
    t0=time.time()
    main()
    print ("Total time elapsed: %f" %(time.time()-t0))
    print ("Maximum memory usage: %f MB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000))
    
        














        





