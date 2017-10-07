import math
import numpy as np


def calc_error(u, uh, x, xh):

    s = 0
    gc = [-math.sqrt(3./5.), 0, math.sqrt(3./5.)]
    w = [5./9., 8./9., 5./9.]

    for i in xrange(len(uh)):
        for j in xrange(3):
            x_loc = (gc[j] + 1.)*(xh[i].max() - xh[i].min())/2. + xh[i].min()

            x_idx = (np.abs(x - x_loc)).argmin()
            u_val = u[0][x_idx]

            xh_idx = (np.abs(xh[i] - x_loc)).argmin()
            uh_val = uh[i][xh_idx]

            s += (math.pow((u_val - uh_val), 2)*((xh[i].max() - xh[i].min())/2.)*w[j])

    return math.sqrt(s)
