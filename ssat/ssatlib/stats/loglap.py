import numpy as np

from scipy.optimize import curve_fit

from ssat.ssatlib import stats

d = 1.0
a = 1.0
b = 1.0


def cdf(x):
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None]
        scalar_input = True

    y = np.zeros(x.size)
    mask = x < d
    y[mask] = a / (a + b) * (x[mask] / d)**b
    mask ^= True
    y[mask] = 1. - b / (a + b) * (d / x[mask])**a

    if scalar_input:
        return np.squeeze(y)
    return y


def fit(xp, fp):
    def f(x, d, a, b):
        g = globals()
        g['d'] = d
        g['a'] = a
        g['b'] = b
        return cdf(x)

    def g(x, d, a, b):
        y = np.zeros((x.size, 3))
        mask = x < d
        c = (x[mask] / d)**b
        y[mask, 0] = -a * b * c / (a * d + b * d)
        y[mask, 1] = b * c / (a + b)**2
        y[mask, 2] = a * c * (-1. + (a + b) * np.log(x[mask] / d)) / (a + b)**2
        mask ^= True
        c = (d / x[mask])**a
        y[mask, 0] = -a * b * c / (a * d + b * d)
        y[mask, 1] = -b * c * (-1. +
                               (a + b) * np.log(d / x[mask])) / (a + b)**2
        y[mask, 2] = -a * c / (a + b)**2
        return y

    popt, pcov = curve_fit(f, xp, fp, p0=np.ones(3), jac=g)
    return stats.rmse(fp, cdf(xp))


def pdf(x):
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None]
        scalar_input = True

    y = np.zeros(x.size)
    c = (a * b) / (d * (a + b))
    mask = x < d
    y[mask] = c * (x[mask] / d)**(b - 1.)
    mask ^= True
    y[mask] = c * (d / x[mask])**(a + 1.)

    if scalar_input:
        return np.squeeze(y)
    return y


def ppf(x):
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[None]
        scalar_input = True

    y = np.zeros(x.size)
    mask = x < a / (a + b)
    y[mask] = d * ((a + b) * x[mask] / a)**(1. / b)
    mask ^= True
    y[mask] = d * ((a + b) * (1. - x[mask]) / b)**(-1. / a)

    if scalar_input:
        return np.squeeze(y)
    return y
