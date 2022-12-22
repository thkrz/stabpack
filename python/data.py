# shape
# mask
# cohesion
# friction angle
# unit weight
# pore pressure
# -> theta_i, theta_r, theta_s, k_sat, beta(2) [van genuchten paramter]
import matplotlib.pyplot as plt
import numpy as np


def mask(xp, fp, grid_size=1):
    assert np.all(np.diff(xp) > 0)
    X = np.arange(xp[0], xp[-1], grid_size)
    m, n = X.size, (fp.max() - fp.min()) // grid_size
    M = np.zeros((m, n), dtype=bool)
    N = int(np.interp(X, xp, fp) - fp.min()) // grid_size
    for i in range(m):
        M[i, : N[i]] = 1
    return M, N
