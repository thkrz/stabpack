import h5py
import numpy as np

from ssat import *
from ssat.ssatlib.hydro import WaterRetentionCurve, conductivity
from ssat.ssatlib.stats import loglap

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
    wrc = WaterRetentionCurve()
    wrc.fit(xp, fp)
    return (res, wrc)


def main():
    with h5py.File(cfg['DATAFILE'], 'r+') as f:
        for name in f['LAYER']:
            g = f['LAYER'][name]
            dset = g['GRANULOMETRY']
            rho = dset.attrs['PARTICLE_DENSITY']
            phi = dset.attrs['POROSITY']
            x = dset[:, 0]
            y = dset[:, 1]
            r, wrc = nimmo(x, y, rho, phi)

            if 'EST' in g:
                del g['EST']
            est = g.create_group('EST')
            est.attrs['A'] = wrc.a
            est.attrs['N'] = wrc.n
            est.attrs['RESIDUAL_WATER'] = r

            est.attrs['HYDRAULIC_CONDUCTIVITY'] = Ks(dset[:, :], loglap.ppf)
    return 0
