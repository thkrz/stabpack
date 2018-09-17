import numpy as np

from ssat import *
from ssat.hydrology import HydraulicConductivity
from ssat.stats import loglap

ALPHA = 1.38
CC = 130.


def ptf(h, rho, phi):
    h = np.asarray(h)
    scalar_input = False
    if h.ndim == 0:
        h = h[None]
        scalar_input = True

    r = np.zeros(h.size)
    e = phi / (1. - phi)
    xi = (2. * e) / (3. * CC**2) * (
        (loglap.a * loglap.b) / (loglap.a + loglap.b) *
        (1. / (2. * np.pi * rho)))**(1. - ALPHA)
    mask = h < 1. / (np.sqrt(xi) * loglap.d**ALPHA)
    r[mask] = (loglap.d**(loglap.b * (ALPHA - 1.)) * xi * h[mask]**2)**(
        1. / (loglap.b * (ALPHA - 1.) - 2. * ALPHA))
    mask ^= True
    r[mask] = (loglap.d**(loglap.a * (1. - ALPHA)) * xi * h[mask]**2)**(
        1. / (loglap.a - (2. + loglap.a) * ALPHA))
    t = phi * loglap.cdf(r)

    if scalar_input:
        return np.squeeze(t)
    return t


def nimmo(x, y, rho, phi):
    loglap.fit(x, y)
    res = ptf(1e6, rho, phi)
    xp = np.logspace(0, 6)
    fp = (ptf(x, rho, phi) - res) / (phi - res)
    hc = HydraulicConductivity()
    hc.fit(xp, fp)


def main():
    with open(cfg['ESTFILE']) as f:
        x = np.fromstring(f.readline(), sep=' ')
        n = x.size
        for line in f:
            y = np.fromstring(line, sep=' ')
            nimmo(x, y[:n], y[n], y[n + 1])
    return False
