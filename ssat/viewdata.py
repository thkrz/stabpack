"""View data"""
import matplotlib.pyplot as plt
import numpy as np

from ssat import *
from ssat.ssatlib.data import Slope


def main():
    slope = Slope(cfg['DATAFILE'], cr=cfg['CRR'])

    fig, ax = plt.subplots()

    ax.plot(slope.omega.x, slope.omega.y, color='black')
    ax.set_axisbelow(True)
    ax.set_aspect(1)
    ax.set_xlim(0, slope.span)
    # ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Length')
    ax.set_ylabel('Height')
    top = slope.omega
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
