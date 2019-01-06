"""Calculate slip surface"""
import numpy as np
import os

from queue import Queue
from threading import Thread

from ssat import *
from ssat.ssatlib.data import Slope
from ssat.ssatlib.stabpack import bezier

result = Queue()


def __calc(k, slope):
    result.put('Test')


def initial_values(s):
    off = int(np.round(cfg['MINLEN'] / cfg['PARADX']))
    if off == 0:
        return None
    x = np.arange(0, s.span, step=cfg['PARADX'])
    y = s.top(x)
    n = x.size
    v = []
    for i in range(n):
        for j in range(i + off, n):
            v.append(bezier.asarc3((x[i], y[i]), (x[j], y[j])))
            v.append(bezier.asline((x[i], y[i]), (x[j], y[j])))
    return np.asarray(v)


def calc(chunk, slope):
    for k in chunk:
        __calc(k, slope)


def main():
    slope = Slope(cfg['DATAFILE'], dim=cfg['DIMENSION'])
    ival = initial_values(slope)
    if ival is None:
        return -1
    num_threads = len(os.sched_getaffinity(0))
    threads = []
    chunks = np.array_split(ival, num_threads)
    for i in range(num_threads):
        t = Thread(target=calc, args=(chunks[i], slope))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
