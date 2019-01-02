"""Calculate slip surface"""
import numpy as np
import os

from queue import Queue
from threading import Thread

from ssat import *
from ssat.ssatlib.data import Slope
from ssat.ssatlib.stabpack import bezier

result = Queue()
# slope = Slope.load(cfg['DATAFILE'])
step = int(np.round(cfg['MINLEN'] / cfg['PARADX']))


def initial_values():
    a = np.array([1, 2])
    b = np.array([10, 8])
    p = bezier.arc(6, a, b, 2.0)
    print(p)


def calc(ctrl_points):
    result.put('Test')


def main():
    if step == 0:
        return -1
    print(bezier)
    ctrl_points = initial_values()
    return 0

    num_threads = len(os.sched_getaffinity(0))
    threads = []
    payload = np.array_split(ctrl_points, num_threads)
    for i in range(num_threads):
        t = Thread(target=calc, args=(payload[i], ))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
