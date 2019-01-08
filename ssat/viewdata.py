"""View data"""
import matplotlib.pyplot as plt
import numpy as np

from ssat import *
from ssat.ssatlib.data import Slope

def main():
    slope = Slope(cfg['DATAFILE'], dim=cfg['DIMENSION'])
    o = slope.stratum(0)
    plt.plot(slope.omega[:, 0], slope.omega[:, 1])
    plt.plot(o.omega[:, 0], o.omega[:, 1])
    plt.show()
    return 0
