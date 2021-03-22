from .parser import parser

def parse_TCRex(filename, out_format='CDR3', sep='\t'):
    """
    Import a dataset following the TCRex input format.
    For mor info about the TCRex format, visit the TCRex documentation:
        https://tcrex.biodatamining.be/instructions/
    
    Parse and extract relevant data columns. These include V gene, J gene and CDR3 amino acid sequence:
        - TRBV_gene
        - TRBJ_gene
        - CDR3_beta
    
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
        Parsed TCRex file.
    """
    
    # Column IDs
    v_col = 'TRBV_gene'
    cdr3_col = 'CDR3_beta'
    j_col = 'TRBJ_gene'
    
    return parser(filename=filename, 
                  v_col=v_col, 
                  cdr3_col=cdr3_col, 
                  j_col=j_col,
                  out_format=out_format,
                  sep=sep)