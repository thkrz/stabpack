import numpy as np

from ssat.ssatlib.stats import loglap


def d(i):
    base = 0.002 if i % 2 == 0 else 0.0063
    return base * 10**i


def phi(d):
    return -np.log2(d)
