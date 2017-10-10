import math
import numpy as np


def calc_error(he, dh, basis, equation, g=0, h=0):

    s = 0
    gc = [-math.sqrt(3./5.), 0., math.sqrt(3./5.)]
    w = [5./9., 8./9., 5./9.]

    x = 0

    for i in xrange(len(he)):
        for j in xrange(3):

            x_loc = (gc[j] + 1.)*(he[i])/2. + x

            u_val = equation(x_loc)

            n2 = basis(2, gc[j])
            n1 = basis(1, gc[j])

            if i < len(he)-1:
                uh_val = n1*dh.item(i) + n2*dh.item(i+1)
            else:
                uh_val = n1*dh.item(i) + n2*g

            dxdc = he[i]/2.

            s += (math.pow((u_val - uh_val), 2.)*dxdc*w[j])

        x += he[i]

    return math.sqrt(s)
