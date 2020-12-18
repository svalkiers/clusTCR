from pandas import read_csv, concat, Series
import os
from typing import Union
import re

MATCHER = re.compile('[A-Z]*')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')


def path_in_data(filename):
    return os.path.join(DATA, filename)


def process_csv(file: Union[str, list],
                cdr3_column_name,
                epitope_column_name=None,
                sep=',',
                cdr3_mapping_func=None):
    """
    Processes a file or list of files and returns the cdr3 and epitope data
    """
    # Make list
    files = [file] if isinstance(file, str) else file
    data = files_to_df(files, sep)
    cdr3 = data[cdr3_column_name]

    if cdr3_mapping_func is not None:
        cdr3 = cdr3.map(cdr3_mapping_func)
    cdr3 = filter_cdr3(cdr3)
    if epitope_column_name is None:
        return cdr3
    else:
        return cdr3, data[epitope_column_name]


def filter_cdr3(data: Series):
    for id in data.index:
        if not MATCHER.fullmatch(data[id]):
            data[id] = None
    return data.dropna().drop_duplicates()


def combine_data(values1, values2):
    values2.index += len(values1)
    return concat([values1, values2])


def files_to_df(files, sep):
    data = None
    for f in files:
        new_data = read_csv(f, sep=sep)
        data = combine_data(data, new_data) if data is not None else new_data
    return data


