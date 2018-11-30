import numpy as np


def rmse(y, yp):
    r = np.sqrt(np.square(yp - y).sum() / y.size)
    rr = r / (np.amax(y) - np.amin(y))
    return (rr, r)
