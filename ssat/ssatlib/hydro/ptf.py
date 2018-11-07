import numpy as np


class Pedotransfer:
    g = 9.81
    rho = 997.
    mu = 8.891e-4
    v = 8.917e-7
    f0 = rho * g / mu
    f1 = g / v

    def __init__(self, n=None, dist=None, ppf=None):
        mask = [n, dist, ppf] == None
        self.n = n
        if dist:
            self.d = dist[:, 0]
            self.w = dist[:, 1]
        self.ppf = ppf

    def phi(self, f):
        return -np.log2(self.ppf(f))


def alyamani(ppf) -> float:
    d50 = ppf(.5)
    d10 = ppf(.1)
    I = d10 - d50 * .2
    K = 1300. * (I + .025 * (d50 - d10))**2
    return K / 86400.


def barr(n, d, w, Cs=1.2) -> float:
    if Cs <= 1. or Cs >= 1.35:
        raise ValueError('Cs must be between 1.0 and 1.35')
    S0 = np.multiply(3. / d, w).sum()
    S = Cs * S0 * (1. - n)
    m = n / S
    K = .2 * f0 * n * m**2
    return K


def beyer(ppf) -> float:
    d60 = ppf(.6)
    d10 = ppf(.1)
    C = d60 / d10
    K = 6e-4 * f1 * np.log(500. / C) * d10**2
    return K


def chapuis(n, ppf) -> float:
    e = n / (1. - n)
    K = 2.4622 * ((d10**2 * e**3) / (1. + e))**0.7825
    return K / 100.


def fair() -> None:
    pass


def harleman(ppf) -> float:
    d10 = ppf(.1)
    K = f0 * d10**2
    return K


def hazen(ppf, C=125.) -> float:
    if C <= 100 or C >= 150:
        raise ValueError('C must be between 100 and 150')
    d10 = ppf(.1) / 10.
    K = C * d10**2
    return K / 100.


def kozeny(n, ppf) -> float:
    d10 = ppf(.1)
    K = 8.3e-4 * f1 * n**3 / (1. - n)**2 * d10**2
    return K


def kozeneycarman(n, ppf) -> float:
    d10 = ppf(.1)
    K = 1. / 180. * f0 * n**3 / (1. - n)**2 * d10**2
    return K


def kruger(n, d, w) -> float:
    de = np.copy(d) * .5
    de[1:] = d[:r.size - 2] + d[1:]
    de = np.divide(de, w)
    K = 4.35e-5 * f1 * n / (1. - n)**2 * de**2
    return K


def navfac(n, ppf) -> float:
    d10 = ppf(.1)
    e = n / (1. - n)
    K = 10.**(1.291 * e - 0.6435) * d10**(10.**(0.5504 - 0.2937 * e))
    return K


def pavchich(ppf) -> float:
    K = .35 * f1 * ppf(.17)**2
    return K


def sauerbrei(n, ppf) -> float:
    K = 3.75e-3 * f1 * n**3 / (1. - n)**2 * ppf(.17)**2
    return K


def slichter(n, ppf) -> float:
    K = .01 * f1 * n**(3.287) * ppf(.1)**2
    return K


def terzaghi(n, ppf, grains='smooth') -> float:
    if grains == 'smooth':
        beta = 10.7e-3
    elif grains == 'coarse':
        beta = 6.1e-3
    else:
        raise ValueError('invalid grains')
    K = beta * f1 * ((n - .13) / np.cbrt(1. - n))**2 * ppf(.1)**2
    return K


def vukovic(ppf) -> float:
    K = 4.8e-4 * f1 * ppf(.2)**2.3
    return K


def zamarin(n, d, w) -> float:
    dg = d[1:]
    dd = d[:d.size - 2]
    de = 3. / 2. * w[0] / d[0] + np.multiply(
        w[1:], np.divide(np.log(dg / dd), dg - dd)).sum()
    de = 1. / de
    K = 8.2e-3 * f1 * n**3 / (1. - n)**2 * de**2
    return K


def zunker(n, d, w, grains='smooth') -> float:
    if grains == 'smooth':
        beta = 2.4e-3
    elif grains == 'coarse':
        beta = 1.4e-3
    elif grains == 'nonuniform':
        beta = 1.2e-3
    elif grains == 'irregular':
        beta = 0.7e-3
    else:
        raise ValueError('invalid grains')

    dg = d[1:]
    dd = d[:d.size - 2]
    de = 3. / 2. * w[0] / d[0] + np.multiply(
        w[1:], np.divide(dg - dd, dg * dd * np.log(dg / dd))).sum()
    de = 1. / de
    K = beta * f1 * n / (1. - n) * de**2


def temperature(deg) -> None:
    global rho, mu, v, f0, f1

    xp = np.array(
        [0.01, 10.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0])
    fp = np.array([(1.7914e-3, 1.7918e-5,
                    999.85), (1.3060e-3, 1.3065e-5,
                              999.70), (1.0016e-3, 1.0035e-5, 998.21),
                   (8.9000e-4, 8.9270e-6,
                    997.05), (7.9720e-4, 8.0070e-6,
                              995.65), (6.5270e-4, 6.5790e-6,
                                        992.22), (5.4650e-4, 5.5310e-6,
                                                  988.04), (4.6600e-4,
                                                            4.7400e-6, 983.20),
                   (4.0350e-4, 4.1270e-6,
                    977.76), (3.5400e-4, 3.6430e-6,
                              971.79), (3.1420e-4, 3.2550e-6, 965.31)])
    rho = np.interp(deg, xp, fp[:, 2])
    mu = np.interp(deg, xp, fp[:, 0])
    v = np.interp(deg, xp, fp[:, 1])
    f0 = rho * g / mu
    f1 = g / v
    return K
