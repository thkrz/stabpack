import numpy as np

from scipy.special import binom


def arc(start, end, depth=None, degree=3):
    d = np.linalg.norm(start - end)
    c = .5 * (start + end)
    print(c)
    exit(0)

    # kappa = 4. / 3. * (np.sqrt(2.) - 1.)
    kappa = .5522847498
    if depth:
        R = depth * .5 + d**2 / (8. * depth)
        phi = np.arccos(1. - depth / R)
    else:
        R = .5 * d
        phi = .25 * np.pi
    m = int(.5 * (degree + 1))
    k = np.zeros((degree + 1, 2))
    k[:m, :] = start
    k[m:, :] = end
    sink = kappa * R * np.sin(phi)
    cosk = kappa * R * np.cos(phi)
    k[1:m, 0] += cosk
    k[1:m, 1] -= sink
    k[m:degree, 0] += cosk
    k[m:degree, 1] -= sink
    print(k)
    return k


def B(t, n, i):
    return binom(n, i) * np.power(t, i) * np.power(1. - t, n - i)


def curve(t, k):
    t = np.asarray(t)
    scalar_input = False
    if t.ndim == 0:
        t = t[None]
        scalar_input = True

    n = k.shape[0] - 1
    b = np.zeros((t.size, k.shape[1]))
    for i in range(n + 1):
        b += np.outer(B(t, n, i), k[i])

    if scalar_input:
        return np.squeeze(b)
    return b


def line(start, end, degree=3):
    n = degree + 1
    t = np.linspace(0, 1, num=n)
    return start + np.outer(t, end - start)


def surface(u, v, k):
    n, m = k.shape
