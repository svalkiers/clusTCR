from .parser import parser

def parse_AIRR(filename, out_format='CDR3', sep='\t'):
    """
    Import a dataset following the AIRR standards data format (MiAIRR Standard Data Representation).
    For mor info about the AIRR format, visit the AIRR Community documentation page:
        https://docs.airr-community.org/en/stable/datarep/format.html
    
    Parse and extract relevant data columns. These include V gene, J gene and CDR3 amino acid sequence:
        - v_call
        - j_call
        - junction_aa
    
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
        Parsed AIRR file.
    """
    
    # Column IDs
    v_col = 'v_call'
    cdr3_col = 'junction_aa'
    j_col = 'j_call'
    
    return parser(filename=filename, 
                  v_col=v_col, 
                  cdr3_col=cdr3_col, 
                  j_col=j_col,
                  out_format=out_format,
                  sep=sep)