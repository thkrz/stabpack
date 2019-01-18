"""View data"""
import matplotlib.pyplot as plt
import numpy as np

from ssat import *
from ssat.ssatlib.data import Slope


def main():
    slope = Slope(cfg['DATAFILE'], cfg['CRR'])

    fig, ax = plt.subplots()

    ax.plot(slope.ridge.x, slope.ridge.y, color='black')
    ax.set_axisbelow(True)
    ax.set_aspect(1)
    ax.set_xlim(*slope.span)
    # ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Length')
    ax.set_ylabel('Height')
    top = slope.ridge
    for stratum in slope:
        if stratum.isdebris():
            ax.fill_between(
                top.x,
                top.y,
                stratum.bottom.y,
                facecolor='none',
                hatch='o',
                edgecolor='black')
            top = stratum.bottom

    plt.grid(True)
    plt.show()
    return 0
