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

result = Queue()


def __calc(k, iscrack, s):
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
        if not all(0 < y - v < maxdepth) or any(u[1:] - u[:u.size - 1] < 0):
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
        return np.minimum(maxmu, mos.razdolsky(slices, F=F))

    p0 = k[1:n, :]
    mu = minimize(f, np.reshape(p0, p0.size), method='Nelder-Mead')
    if mu < 1.0:
        result.put((mu, k))


def initial_values(s, num):
    off = int(np.round(cfg['MINLEN'] / cfg['PARADX']))
    if off == 0:
        return None
    x = np.arange(0, s.span, step=cfg['PARADX'])
    y = s.top(x)
    n = x.size
    v = []
    for i in range(n):
        iscrack = s.iscrack(x[i])
        lyr = s.stratum((x[i], y[i]))
        y[i] -= lyr.depth()
        for j in range(i + off, n):
            alpha = np.arctan((y[j] - y[i]) / (x[j] - x[i]))
            if np.degrees(alpha) < cfg['MINDESC']:
                continue
            k = bezier.asarc3((x[i], y[i]), (x[j], y[j]))
            v.append((k, iscrack))
            k = bezier.asline((x[i], y[i]), (x[j], y[j]))
            v.append((k, iscrack))
    return (v[i:i + num] for i in range(0, len(v), num))


def calc(chunk, s):
    for k, iscrack in chunk:
        __calc(k, iscrack, s)


def main():
    num_threads = len(os.sched_getaffinity(0))
    slope = Slope(cfg['DATAFILE'], dim=cfg['DIMENSION'])
    chunks = initial_values(slope, num_threads)
    if chunks is None:
        return -1
    threads = []
    for i in range(num_threads):
        t = Thread(target=calc, args=(chunks[i], slope))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
