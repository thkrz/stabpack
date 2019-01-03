import numpy as np

from scipy.special import binom


def asarc(start, end):
    start = np.asarray(start)
    end = np.asarray(end)

    dx = end - start
    theta = .5 * np.pi - np.arctan(dx[1] / dx[0])
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    a = start
    b = np.dot(rot, dx) + a

    kappa = .5522847498
    phi = .25 * np.pi
    R = np.abs(a[1] - b[1]) / (2. * np.sin(phi))
    k = np.zeros((4, 2))
    k[0, :] = a
    k[1, 0] = a[0] + kappa * R * np.sin(phi)
    k[1, 1] = a[1] + kappa * R * np.cos(phi)
    k[2, 0] = b[0] + kappa * R * np.sin(phi)
    k[2, 1] = b[1] - kappa * R * np.cos(phi)
    k[3, :] = b
    rrot = np.linalg.inv(rot)
    for i in range(1, 4):
        k[i, :] = np.dot(rrot, k[i, :] - a) + a
    return k


def asline(start, end, degree=3):
    start = np.asarray(start)
    end = np.asarray(end)
    t = np.linspace(0, 1, num=(degree + 1))
    return start + np.outer(t, end - start)


def bernstein(t, n, i):
    return binom(n, i) * np.power(t, i) * np.power(1. - t, n - i)


def curve(k, t):
    t = np.asarray(t)
    scalar_input = False
    if t.ndim == 0:
        t = t[None]
        scalar_input = True

    n = k.shape[0] - 1
    b = np.zeros((t.size, k.shape[1]))
    for i in range(n + 1):
        b += np.outer(bernstein(t, n, i), k[i])

    if scalar_input:
        return np.squeeze(b)
    return b
