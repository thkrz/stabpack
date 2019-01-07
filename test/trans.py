"""Derive parameters"""
import h5py
import matplotlib.pyplot as plt
import numpy as np

from ssat import *
from ssat.ssatlib.constants import G
from ssat.ssatlib import hydro
from ssat.ssatlib.pedo import diaspace, grade, sort
from ssat.ssatlib.stats import loglap

phi = 0.5
rho = 2650.
d = None
w = None

F0 = hydro.rho() * G / hydro.mu()
F1 = F0


def ks_alyamani():
    d50 = loglap.ppf(.5)
    d10 = loglap.ppf(.1)
    I = d10 - d50 * .2
    K = 1300. * (I + .025 * (d50 - d10))**2
    return K * 1.1574074e-05


def ks_barr(Cs=1.2):
    if Cs <= 1. or Cs >= 1.35:
        raise ValueError('Cs must be between 1.0 and 1.35')
    S0 = np.multiply(3. / (d * .001), w).sum()
    S = Cs * S0 * (1. - phi)
    m = phi / S
    K = .2 * F0 * phi * m**2
    return K


def ks_beyer():
    d60 = loglap.ppf(.6)
    d10 = loglap.ppf(.1)
    if d10 <= 0.06 or d10 >= 0.6:
        return None
    d60 *= 0.001
    d10 *= 0.001
    C = d60 / d10
    if C <= 1. or C >= 20.:
        return None
    K = F0 * 6e-4 * np.log(500. / C) * d10**2
    return K


def ks_chapuis():
    d10 = loglap.ppf(.1)
    if d10 <= 0.03 or d10 >= 3.0:
        return None
    e = phi / (1. - phi)
    if e <= 0.3 or e >= 0.7:
        return None
    K = 2.4622 * ((d10**2 * e**3) / (1. + e))**0.7825
    return K * .01


# def ks_fair() -> None:
#     pass


def ks_harleman():
    d10 = loglap.ppf(.1) * 0.001
    K = F0 * d10**2
    return K


def ks_hazen(C=125.):
    if C <= 100 or C >= 150:
        raise ValueError('C must be between 100 and 150')
    d10 = loglap.ppf(.1) * .1
    K = C * d10**2
    return K * .01


def ks_kozeny():
    d10 = loglap.ppf(.1) * 0.001
    K = 8.3e-4 * F1 * phi**3 / (1. - phi)**2 * d10**2
    return K


def ks_kozeneycarman():
    d10 = loglap.ppf(.1)
    if d10 >= 3.0:
        return None
    d10 *= 0.001
    K = 1. / 180. * F0 * phi**3 / (1. - phi)**2 * d10**2
    return K


def ks_kruger():
    de = np.copy(d) * .5
    de[1:] = d[:d.size - 1] + d[1:]
    de = 1. / np.divide(w, de * 0.001).sum()
    K = 4.35e-5 * F1 * phi / (1. - phi)**2 * de**2
    return K


def ks_navfac():
    d10 = loglap.ppf(.1)
    if d10 <= 0.1 or d10 >= 2.0:
        return None
    if d10 / loglap.ppf(.05) < 1.4:
        return None
    e = phi / (1. - phi)
    if e <= 0.3 or e >= 0.7:
        return None
    K = 10.**(1.291 * e - 0.6435) * d10**(0.5504 - 0.2937 * e)
    return K


def ks_pavchich():
    d17 = loglap.ppf(.17)
    if d17 <= 0.06 or d17 >= 1.5:
        return None
    d17 *= 0.001
    K = .35 * F1 * d17**2
    return K


def ks_sauerbrei():
    d17 = loglap.ppf(.17)
    if d17 >= 0.5:
        return None
    d17 *= 0.001
    K = 3.75e-3 * F1 * phi**3 / (1. - phi)**2 * d17**2
    return K


def ks_slichter():
    d10 = loglap.ppf(.1)
    if d10 <= 0.01 or d10 >= 5.0:
        return None
    d10 *= 0.001
    K = .01 * F1 * phi**(3.287) * d10**2
    return K


def ks_terzaghi(grains='smooth'):
    if grains == 'smooth':
        beta = 10.7e-3
    elif grains == 'coarse':
        beta = 6.1e-3
    else:
        raise ValueError('invalid grains')
    d10 = loglap.ppf(.1) * 0.001
    K = beta * F1 * ((phi - .13) / np.cbrt(1. - phi))**2 * d10**2
    return K


def ks_vukovic():
    d20 = loglap.ppf(.2) * 0.001
    K = 4.8e-4 * F1 * d20**2.3
    return K


def ks_zamarin():
    dg = d[1:]
    dd = d[:d.size - 1]
    de = 3. / 2. * w[0] / d[0] + np.multiply(
        w[1:], np.divide(np.log(dg / dd), dg - dd)).sum()
    de = 1. / de
    de *= 0.001
    K = 8.2e-3 * F1 * phi**3 / (1. - phi)**2 * de**2
    return K


def ks_zunker(grains='smooth'):
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
    dd = d[:d.size - 1]
    de = 3. / 2. * w[0] / d[0] + np.multiply(
        w[1:], np.divide(dg - dd, dg * dd * np.log(dg / dd))).sum()
    de = 1. / de
    de *= 0.001
    K = beta * F1 * phi / (1. - phi) * de**2


def t_nimmo(h, alpha=1.38, cc=130.):
    h = np.asarray(h)
    scalar_input = False
    if h.ndim == 0:
        h = h[None]
        scalar_input = True

    r = np.zeros(h.size)
    e = phi / (1. - phi)
    xi = (2. * e) / (3. * cc**2) * (
        (loglap.a * loglap.b) / (loglap.a + loglap.b) *
        (1. / (2. * np.pi * rho)))**(1. - alpha)
    mask = h < 1. / (np.sqrt(xi) * loglap.d**alpha)
    r[mask] = (loglap.d**(loglap.b * (alpha - 1.)) * xi * h[mask]**2)**(
        1. / (loglap.b * (alpha - 1.) - 2. * alpha))
    mask ^= True
    r[mask] = (loglap.d**(loglap.a * (1. - alpha)) * xi * h[mask]**2)**(
        1. / (loglap.a - (2. + loglap.a) * alpha))
    t = phi * loglap.cdf(r)

    if scalar_input:
        return np.squeeze(t)
    return t


# def fitshit(x, y, rho, phi):
#     loglap.fit(x, y)
#     res = ptf(1e6, rho, phi)
#     xp = np.logspace(0, 6)
#     fp = (ptf(x, rho, phi) - res) / (phi - res)
#     wrc = WaterRetentionCurve()
#     wrc.fit(xp, fp)
#     return (res, wrc)


def plot():
    fig, ax = plt.subplots()
    ax.set_xscale('log', basex=2)

    ax.scatter(d, np.cumsum(w))
    x = np.linspace(0.0063, 2, num=100)
    ax.plot(x, loglap.cdf(x))
    plt.show()


def main():
    global phi, rho, d, w
    ks = []
    header = 'grade;sort;'
    for k, v in globals().items():
        if k.startswith('ks_'):
            header += k[3:] + ';'
            ks.append(v)
    arr = np.loadtxt(sys.argv[1])
    d = diaspace()
    print(header[:-1])
    for p in arr:
        rho_b = p[0]
        rho = p[1]
        w = p[2:]
        if w[0] < 0.001:
            print("-;")
            continue
        rrmse, rmse = loglap.fit(d, np.cumsum(w))
        log('rrms/rms: {:.2f}/{:.2f}'.format(rrmse, rmse))
        if rho_b < 0.1:
            phi = .255 * (1. + .83**(loglap.ppf(.6) / loglap.ppf(.1)))
        else:
            phi = 1. - (rho_b / rho)
        res = []
        for fcn in ks:
            res.append(str(fcn() or -999.))
        print(';'.join([str(phi), grade(w), sort(loglap.ppf), *res]))
    return 1

    arr = np.loadtxt(sys.argv[1])
    d = arr[::-1, 0]
    w = arr[::-1, 1] / 100.
    w[1:] -= w[:-1].copy()
    rrmse, rmse = loglap.fit(d, np.cumsum(w))
    log('rrms/rms: {:.2f}/{:.2f}'.format(rrmse, rmse))
    log('d10: {:.4f}'.format(float(loglap.ppf(.1))))
    log('d60: {:.4f}'.format(float(loglap.ppf(.6))))
    # n_mess = 1. - (rho_b / rho_p)
    n_calc = .255 * (1. + .83**(loglap.ppf(.6) / loglap.ppf(.1)))
    log('n: {:.1f}%'.format(n_calc * 100.))
    phi = n_calc
    log('grade: ' + grade(w))
    log('sorting: ' + sort(loglap.ppf))
    log('Estimation of Hydraulic Conductivity:')
    for model, fcn in ks.items():
        print('{:14s}: {:.4e} m/s'.format(model[3:], fcn() or -999.))
    # with h5py.File(cfg['DATAFILE'], 'r') as f:
    #     for name in f['LAYER']:
    #         g = f['LAYER'][name]
    #         # rho = g.attrs['PARTICLE_DENSITY']
    #         # phi = g.attrs['POROSITY']
    #         dset = g['GRANULOMETRY']
    #         d = dset[:, 0]
    #         w = dset[:, 1]
    #         log('Fitting grain size distribution')
    #         rrmse, rmse = loglap.fit(d, np.cumsum(w))
    #         phi = .255 * (1. + .83**(loglap.ppf(.6) / loglap.ppf(.1)))
    #         log('rrms/rms: {:.4f}/{:.4f}'.format(rrmse, rmse))
    #         log('grading: ' + grade(w))
    #         log('sorting: ' + sort(loglap.ppf))
    #         res = []
    #         for fcn in ks.values():
    #             res.append(str(fcn() or -999.))
    #         print(';'.join([grade(w), name, *res]))

    # r, wrc = nimmo(x, y, rho, phi)
    return 1
