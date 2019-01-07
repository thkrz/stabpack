import numpy as np

from collections import namedtuple

Slice = namedtuple('Slice', ['alpha', 'b', 'c', 'gamma', 'h', 'phi', 'u'])


def razdolsky(slices, crest=0.0, F=1.3):
    n = len(slices) + 1
    E = np.zeros(n)
    T = np.zeros(n)
    E[0] = crest
    T[0] = crest
    for i, S in enumerate(slices):
        tana = np.tan(S.alpha)
        tanp = np.tan(S.phi) / F
        tan2 = tana * tana
        c = S.c / F
        E[i + 1] = 2. * T[i] * (tana - tanp) + E[i] * (
            1. + 2. * tana * tanp - tan2) + S.gamma * S.b / 3. * (
                S.h[1] + 2. * S.h[0]) * (tana - tanp) - S.b * (1. + tan2) * (
                    c - S.u * tanp)
        E[i + 1] /= 1. + tan2
        T[i + 1] = (E[i + 1] + E[i]) * tana - T[i] + S.gamma * S.b / 6. * (
            S.h[1] - S.h[0])
    return 1. - E[n - 1] / np.amax(E)
