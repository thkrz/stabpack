import numpy as np

from scipy.optimize import curve_fit

from ssat.ssatlib import stats
from ssat.ssatlib.constants import G

temperature = 20.


class Retention:
    def __init__(self, a=0.0005, n=2.0):
        self.a = a
        self.m = 1. - 1. / n
        self.n = n

    def fit(self, xp, fp):
        def f(x, a, n):
            self.__init__(a, n)
            return self.Theta(x)

        def g(x, a, n):
            y = np.zeros((x.size, 2))
            m = 1. - 1. / n
            p = (a * x)**n
            c = -p * (1. / (1. + p))**(1. + m) * m
            y[:, 0] = c * n / a
            y[:, 1] = c * np.log(a * x)
            return y

        popt, pcov = curve_fit(f, xp, fp, p0=[0.0005, 2.0], jac=g)
        return stats.rmse(fp, self.Theta(xp))

    def psi(self, T):
        T = np.asarray(T)
        scalar_input = False
        if T.ndim == 0:
            T = T[None]
            scalar_input = True

        p = (1. / self.a**self.n * (1. / T)**(1. / self.m) - 1.)**(1. / self.n)
        pe = self.a * (0.046 * self.m + 2.07 * self.m**2 + 19.5 * self.m**3
                       ) / (1. + 4.7 * self.m + 16. * self.m**2)
        p[p < pe] = pe

        if scalar_input:
            return np.squeeze(p)
        return p

    def K_rel(self, T):
        return np.sqrt(T) * (1. - (1. - T**(1. / self.m))**self.m)**2

    def Theta(self, h):
        return (1. / (1. + (self.a * h)**self.n))**self.m


def K(k):
    k * rho() * G / mu()


def mu(T=temperature):
    T += 273.15
    if T < 273 or T > 373:
        raise ValueError
    A, B, C = -3.7188, 578.919, -137.546
    return np.exp(A + B / (C + T))


def rho(T=temperature):
    T += 273.15
    if T < 273 or T > 648:
        raise ValueError
    A, B, C, D = 0.14395, 0.0112, 649.727, 0.05107
    return A / np.power(B, 1. + np.power(1. - T / C, D))
