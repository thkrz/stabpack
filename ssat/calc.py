"""Calculate slip surface"""
import numpy as np
import os

from queue import Queue
from threading import Thread

from ssat import *
from ssat.ssatlib.data import Slope

result = Queue()


def calc(i):
    result.put(i)


def main():
    s = Slope.load(cfg['DATAFILE'])
    num_threads = len(os.sched_getaffinity(0))
    threads = []
    for i in range(num_threads):
        t = Thread(target=calc, args=(i,))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    while not result.empty():
        print(result.get())
    return 0
