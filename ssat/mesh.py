import numpy as np

from ssat import *

def mesh(relief, strata):
    pass


def main():
    with open(cfg['SHAPEFILE']) as f:
        lines = []
        for line in in f:
            line = line.strip()
            if line and not line.startswith('#'):
                lines.append(line)
    dim = 2
    try:
        i = 0
        j = int(lines[i])
        i += 1
        r = np.fromstring('\n'.join(lines[i:i + j], sep=' ')).reshape((dim, -1))
        i = j
        j = int(lines[i])
        i += 1
        s = np.fromstring('\n'.join(lines[i:i + j], sep=' ')).reshape((dim, -1))
    except:
        die('invalid shapefile.')
    mesh(r, s)
    return False
