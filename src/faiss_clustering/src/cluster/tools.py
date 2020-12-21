import numpy as np
import pandas as pd
import faiss

from ..profile.profile import make_profile


def make_profiles(values: pd.Series,
                  properties,
                  vector_mapping_func=None,
                  add_constant=False,
                  extra_feature_size=0,
                  add_average=False,
                  add_length=False):
    matrix, profile_length = make_matrix(values, properties, extra_feature_size, add_constant, add_average, add_length)
    for i, elem in enumerate(values):
        vec = make_vec(elem, profile_length, properties, vector_mapping_func, add_average, add_constant, add_length)
        matrix[i] = vec
    return matrix


def make_matrix(values, properties, extra_feature_size, add_constant, add_average, add_length):
    n = len(values)
    profile_length = values.str.len().max()
    y_size = profile_length * len(properties) + extra_feature_size
    if add_constant:
        y_size += 1
    if add_length:
        y_size += 1
    if add_average:
        y_size += len(properties)
    matrix = np.zeros((n, y_size)).astype('float32')
    return matrix, profile_length


def make_vec(sequence, max_length, properties, vector_mapping_func=None, add_average=False, add_constant=False, add_length=False):
    vec = []
    for prop in properties:
        profile = make_profile(sequence, prop)
        if vector_mapping_func is not None:
            profile = vector_mapping_func(profile)
        if add_average:
            profile.append(sum(profile))
        vec.extend(pad_vector(profile, max_length))
    if add_constant:
        vec.append(0)
    if add_length:
        vec.append(len(sequence))
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


# The following line disables an inspection in PyCharm, the train and add methods only require 1 argument instead of 2
# noinspection PyArgumentList
def cluster_with_faiss(matrix, items_per_cluster, ids=None):
    ncentroids = matrix.shape[0] // items_per_cluster
    if ncentroids == 0:
        ncentroids = 1
    dimension = matrix.shape[1]
    quantizer = faiss.IndexFlatL2(dimension)
    index = faiss.IndexIVFFlat(quantizer, dimension, ncentroids)
    index.train(matrix)
    if ids is not None:
        index.add_with_ids(matrix, ids)
    else:
        index.add(matrix)
    return index

