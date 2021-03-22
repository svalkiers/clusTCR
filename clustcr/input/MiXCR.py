from .parser import parser

def parse_MiXCR(filename, out_format='CDR3', sep='\t'):
    
    """
    Import a dataset in the MiXCR data format.
    For mor info about the MiXCR format, visit the MiXCR documentation page:
        https://mixcr.readthedocs.io/en/master/
    
    Parse and extract relevant data columns. These include V gene, J gene and CDR3 amino acid sequence:
        - allVHitsWithScore
        - allJHitsWithScore
        - aaSeqCDR3
    
    Arguments
    ---------
    filename: str
        Name of the input file.
    out_format: str, default = CDR3
        The output format. Available output formats include:
        CDR3, GLIPH2, TCRDist
    sep: str, default = \t
        Column separator.
    
    Returns
    -------
    data: pd.DataFrame
        Parsed MiXCR file.
    """
    
    v_col = 'allVHitsWithScore'
    cdr3_col = 'aaSeqCDR3'
    j_col = 'allJHitsWithScore'
    
    return parser(filename=filename, 
                  v_col=v_col, 
                  cdr3_col=cdr3_col, 
                  j_col=j_col,
                  out_format=out_format,
                  sep=sep)