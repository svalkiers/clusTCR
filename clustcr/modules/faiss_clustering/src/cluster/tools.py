import numpy as np
import pandas as pd
import parmap
import multiprocessing
from ..profile.profile import make_profile


def make_profiles(values: pd.Series, properties, size, n_cpus):
    matrix, profile_length = make_matrix(values, properties, size)
    with multiprocessing.Pool(n_cpus) as pool:
        vecs = parmap.map(make_vec,
                          values,
                          profile_length,
                          properties,
                          pm_parallel=True,
                          pm_pool=pool)
    for i in range(len(vecs)):
        matrix[i] = vecs[i]
    return matrix


def make_matrix(values, properties, profile_length):
    n = len(values)
    if profile_length is None:
        profile_length = values.str.len().max()
    y_size = profile_length * len(properties)
    matrix = np.zeros((n, y_size)).astype('float32')
    return matrix, profile_length


def make_vec(sequence, max_length, properties):
    vec = []
    for prop in properties:
        profile = make_profile(sequence, prop)
        vec.extend(pad_vector(profile, max_length))
    return vec


def pad_vector(vec, n):
    padding = n - len(vec)
    half = padding // 2
    right = half if padding % 2 == 0 else half + 1
    return [0] * half + vec + [0] * right
