"""Calculate slip surface"""
import numpy as np
import os

from queue import Queue
from threading import Thread

from ssat import *
from ssat.ssatlib.data import Slope
from ssat.ssatlib.stabpack import bezier

result = Queue()
slope = None

def __calc(k):
    result.put('Test')


def initial_values():
    step = int(np.round(cfg['MINLEN'] / cfg['PARADX']))
    if step == 0:
        return None
    x = np.arange(0, 30, step=cfg['PARADX'])
    y = np.ones(x.size)
    n = x.size
    iv = []
    for i in range(n):
        for j in range(i + step, n):
            iv.append(bezier.asarc((x[i], y[i]), (x[j], y[j])))
            iv.append(bezier.asline((x[i], y[i]), (x[j], y[j])))
    return iv


def calc(chunk):
    for k in chunk:
        __calc(k)


def main():
    num_threads = len(os.sched_getaffinity(0))
    threads = []
    chunks = np.array_split(initial_values(), num_threads)
    print(chunks)
    return 0
    for i in range(num_threads):
        t = Thread(target=calc, args=(chunks[i], ))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
