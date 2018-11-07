import numpy as np

from scipy.optimize import curve_fit


class WaterRetentionCurve:
    def __init__(self, a=0.0005, n=2.0):
        self.a = a
        self.m = 1. - 1. / n
        self.n = n

    def fit(self, xp, fp):
        def f(x, a, n):
            self.__init__(a, n)
            return self.theta(x)

        def grad(x, a, n):
            y = np.zeros((x.size, 2))
            m = 1. - 1. / n
            p = (a * x)**n
            c = -p * (1. / (1. + p))**(1. + m) * m
            y[:, 0] = c * n / a
            y[:, 1] = c * np.log(a * x)
            return y

        popt, pcov = curve_fit(f, xp, fp, p0=[0.0005, 2.0], jac=grad)

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

    def theta(self, h):
        return (1. / (1. + (self.a * h)**self.n))**self.m
