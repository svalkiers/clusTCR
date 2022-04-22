""" The chemical properties that can be used for the clustering """
from itertools import chain, combinations

def powerset(s):
    """ Powerset without the empty tuple """
    pset = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    return list(filter(None, pset))


BASICITY = 'basicity'
HELICITY = 'helicity'
HYDROPHOBICITY = 'hydrophobicity'
MUTATION_STABILITY = 'mutation_stability'

ATCHLEY_F1 = 'atchley1'
ATCHLEY_F2 = 'atchley2'
ATCHLEY_F3 = 'atchley3'
ATCHLEY_F4 = 'atchley4'
ATCHLEY_F5 = 'atchley5'
ATCHLEY_FACTORS = [ATCHLEY_F1, ATCHLEY_F2, ATCHLEY_F3, ATCHLEY_F4, ATCHLEY_F5]

Z1 = 'z_scores1'
Z2 = 'z_scores2'
Z3 = 'z_scores3'
Z_SCORES = [Z1, Z2, Z3]

ALL_PROPS = [BASICITY, HELICITY, HYDROPHOBICITY, MUTATION_STABILITY, *ATCHLEY_FACTORS, *Z_SCORES]
OPTIMAL = [MUTATION_STABILITY, *Z_SCORES]
PROPERTY_COMBINATIONS = powerset(ALL_PROPS)
