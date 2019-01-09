"""View data"""
import matplotlib.pyplot as plt
import numpy as np

from ssat import *
from ssat.ssatlib.data import Slope

def main():
    slope = Slope(cfg['DATAFILE'], dim=cfg['DIMENSION'], cr=cfg['CRR'])

    plt.plot(slope.omega.x, slope.omega.y)
    for stratum in slope:
        if stratum.isdebris():
            plt.plot(stratum.omega.x, stratum.omega.y)
    plt.show()
    return 0
