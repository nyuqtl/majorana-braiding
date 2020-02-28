import numpy as np


def are_close(val1, val2, atol=1e-08):
    return np.isclose([val1], [val2], atol=atol)[0]


def osum(lst):
    return np.sum(np.array(lst, dtype=object))
