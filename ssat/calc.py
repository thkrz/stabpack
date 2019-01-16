"""Calculate slip surface"""
import numpy as np
import os

from queue import Queue
from threading import Thread
from scipy.optimize import minimize

from ssat import *
from ssat.ssatlib import hydro
from ssat.ssatlib.constants import G
from ssat.ssatlib.data import Slope
from ssat.ssatlib.stabpack import bezier, mos

import matplotlib.pyplot as plt

result = Queue()


def __calc(k, f0, s):
    dim = cfg['DIMENSION']
    n = cfg['DEGREE']
    F = cfg['FOS']
    maxmu = cfg['MAXMU']
    maxdepth = np.linalg.norm(k[0, :] - k[n, :]) * cfg['MAXDEPTH']
    t = np.linspace(0, 1, num=cfg['NPARAM'])

    # Spaghetti code ^ 3
    # TODO: u and gamma, crest/cracks
    def f(x, *args):
        k[1:n, :] = x.reshape((n - 1, dim))
        u, v = np.hsplit(bezier.curve(k, t), 2)
        y = s.top(u)
        dy = y - v
        if np.logical_or(dy < 0, dy > maxdepth).any() or np.any(u[1:] - u[:u.size - 1] < 0):
            return maxmu
        slices = []
        for i in range(u.size - 1):
            b = u[i + 1] - u[i]
            alpha = np.arctan((v[i + 1] - v[i]) / b)
            mx = .5 * (u[i] + u[i + 1])
            my = .5 * (v[i] + v[i + 1])
            h = np.abs((y[i] - v[i], y[i + 1] - v[i + 1]))
            layer = s.stratum((mx, my))
            c = layer.c
            gamma = s.gamma((mx, my))
            phi = layer.phi
            u = layer.u(my)
            slices.append(mos.Slice(alpha, b, c, gamma, h, phi, u))
        return np.minimum(maxmu, mos.razdolsky(slices, crest=f0, F=F))

    p0 = k[1:n, :]
    mu = minimize(f, np.reshape(p0, p0.size), method='Nelder-Mead')
    if mu < 1.0:
        result.put((mu, k))


def boundary_values(s, num):
    a = cfg['MINLEN']
    beta = cfg['MINDESC']
    dx = cfg['PARADX']
    if a < dx or any([a < 0, beta < 0, dx < 0]):
        raise ValueError

    off = int(np.round(a / dx))
    if off == 0:
        return None
    x = np.arange(*s.span, step=dx)
    y = s.top(x)
    n = x.size
    b = []
    eps = dx * .5
    for i in range(n):
        f0 = 0.0
        if s.iscrack(x[i], eps):
            l = s.stratum((x[i], y[i]))
            h = l.depth()
            print(h)
            y[i] -= h
            f0 = hydro.rho() * G * h
        for j in range(i + off, n):
            alpha = np.arctan((y[i] - y[j]) / (x[j] - x[i]))
            if np.degrees(alpha) < beta:
                continue
            b.append(((x[i], y[i]), (x[j], y[j]), f0))
    return np.array_split(np.asarray(b), num)


def calc(chunk, s):
    for x0, x1, f0 in chunk:
        __calc(bezier.asarc3(x0, x1), f0, s)
        __calc(bezier.asline(x0, x1), f0, s)


def main():
    slope = Slope(cfg['DATAFILE'], cr=cfg['CRR'])
    num_threads = len(os.sched_getaffinity(0))
    chunks = boundary_values(slope, num_threads)
    if chunks is None:
        return -1
    threads = []
    for chunk in chunks:
        t = Thread(target=calc, args=(chunk, slope))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
