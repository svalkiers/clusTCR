import numpy as np
import pandas as pd
import parmap
import multiprocessing
import time
from ..profile.profile import make_profile


def make_profiles(values: pd.Series, properties, size, n_cpus):
    a = time.time()
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
    print(time.time() - a)
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
    add_left = False
    while len(vec) != n:
        if add_left:
            vec.insert(0, 0)
        else:
            vec.append(0)
        add_left = not add_left
    return vec

