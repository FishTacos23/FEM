import math
import numpy as np


def calc_error(uh, xh, dh, basis, equation, p, g=0):

    s = 0
    gc = [-math.sqrt(3./5.), 0., math.sqrt(3./5.)]
    w = [5./9., 8./9., 5./9.]

    for i in xrange(len(uh)):
        for j in xrange(3):

            x_loc = (gc[j] + 1.)*(xh[i].max() - xh[i].min())/2. + xh[i].min()

            u_val = equation(x_loc)

            n2 = basis(2, gc[j])
            n1 = basis(1, gc[j])

            if i < len(uh)-1:
                uh_val = n1*dh.item(i) + n2*dh.item(i+1)
            else:
                uh_val = n1*dh.item(i) + n2*g

            dxdc = (xh[i].max() - xh[i].min())/2.

            s += (math.pow((u_val - uh_val), 2.)*dxdc*w[j])

            if p:
                print str(uh_val) + ', ' + str(u_val)

    if p:
        str(math.sqrt(s))

    return math.sqrt(s)
