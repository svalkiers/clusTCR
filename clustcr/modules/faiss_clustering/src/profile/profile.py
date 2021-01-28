# AUTHOR: Pieter Meysman
# Slightly altered by Max Van Houcke
# - In the make_profile a default value of 0 is set for chemical properties
# - Added mutation stability and pI
# - Added Atchley Factors
# - Added z_scores
# - Removed unused functions

from scipy.stats import zscore
from pyteomics.electrochem import pI
from .atchley import ATCHLEY_FACTORS
from .z_scores import Z_SCORES


# Z-normalization of dictionaries
def znorm_scipy(d):
    keys, vals = zip(*d.items())
    return dict(zip(keys, zscore(vals, ddof=1)))


basicity = {'A': 206.4, 'B': 210.7, 'C': 206.2, 'D': 208.6, 'E': 215.6, 'F': 212.1, 'G': 202.7,
            'H': 223.7, 'I': 210.8, 'K': 221.8, 'L': 209.6, 'M': 213.3, 'N': 212.8, 'P': 214.4,
            'Q': 214.2, 'R': 237.0, 'S': 207.6, 'T': 211.7, 'V': 208.7, 'W': 216.1, 'X': 210.2,
            'Y': 213.1, 'Z': 214.9}

hydrophobicity = {'A': 0.16, 'B': -3.14, 'C': 2.50, 'D': -2.49, 'E': -1.50, 'F': 5.00, 'G': -3.31,
                  'H': -4.63, 'I': 4.41, 'K': -5.00, 'L': 4.76, 'M': 3.23, 'N': -3.79, 'P': -4.92,
                  'Q': -2.76, 'R': -2.77, 'S': -2.85, 'T': -1.08, 'V': 3.02, 'W': 4.88, 'X': 4.59,
                  'Y': 2.00, 'Z': -2.13}

helicity = {'A': 1.24, 'B': 0.92, 'C': 0.79, 'D': 0.89, 'E': 0.85, 'F': 1.26, 'G': 1.15, 'H': 0.97,
            'I': 1.29, 'K': 0.88, 'L': 1.28, 'M': 1.22, 'N': 0.94, 'P': 0.57, 'Q': 0.96, 'R': 0.95,
            'S': 1.00, 'T': 1.09, 'V': 1.27, 'W': 1.07, 'X': 1.29, 'Y': 1.11, 'Z': 0.91}

mutation_stability = {'A': 13, 'C': 52, 'D': 11, 'E': 12, 'F': 32, 'G': 27, 'H': 15, 'I': 10,
                      'K': 24, 'L': 34, 'M':  6, 'N':  6, 'P': 20, 'Q': 10, 'R': 17, 'S': 10,
                      'T': 11, 'V': 17, 'W': 55, 'Y': 31}

elektrochem = {'basicity': znorm_scipy(basicity), 'hydrophobicity': znorm_scipy(hydrophobicity),
               'helicity': znorm_scipy(helicity), 'mutation_stability': znorm_scipy(mutation_stability)}


def make_profile(sequence, prop='basicity'):
    if prop == 'pI':
        return [pI(x) for x in sequence]
    elif 'atchley' in prop:
        # atchley1 for example gives the first factor
        index = int(prop[-1]) - 1
        return [ATCHLEY_FACTORS.get(x, [0, 0, 0, 0, 0])[index] for x in sequence]
    elif 'z_scores' in prop:
        # z_scores1 for example gives the first factor
        index = int(prop[-1]) - 1
        return [Z_SCORES.get(x, [0, 0, 0])[index] for x in sequence]
    else:
        return [elektrochem[prop].get(x, 0) for x in sequence]

