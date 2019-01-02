import numpy as np

from scipy.special import binom


def arc(start, end, depth=None, rank=3):
    kappa = 4. / 3. * (np.sqrt(2.) - 1.)
    d = np.linalg.norm(start - end)
    if depth:
        R = depth * .5 + d**2 / (8. * depth)
        phi = np.arccos(1. - depth / R)
    else:
        R = .5 * d
        phi = .5 * np.pi
    m = int(.5 * (rank + 1))
    k = np.zeros((rank + 1, 2))
    k[:m, :] = start
    k[m:, :] = end
    k[1:rank, 0] += kappa * R * np.sin(phi)
    k[1:rank, 1] -= kappa * R * np.cos(phi)
    return k


def B(t, n, i):
    return binom(n, i) * np.power(t, i) * np.power(1. - t, n - i)


def curve(t, k):
    t = np.asarray(t)
    scalar_input = False
    if t.ndim == 0:
        t = t[None]
        scalar_input = True

    n, m = k.shape
    b = np.zeros((t.size, m))
    for i in range(n):
        b += np.outer(B(t, n, i), k[i])

    if scalar_input:
        return np.squeeze(b)
    return b


def line(start, end, rank=3):
    pass

def surface(u, v, k):
    n, m = k.shape
