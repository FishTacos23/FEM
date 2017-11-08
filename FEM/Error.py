import math
import numpy as np
import Finite_Element_Method as fem


def calc_error(he, dh, basis, equation, model, g=0, h=0, p=2):

    s = 0
    element = model.n
    gc = model.qs
    w = model.ws

    for e in xrange(1, element+1):
        for j, qs in enumerate(gc):

            dnx, ddnx, jac, x, ne = model._basis_x(e, qs)

            u_val = equation(x)

            uh_val = 0

            for a in xrange(1, p + 2):  # loop over the local basis
                uh_val += ne[a-1]*dh[model._ien(e, a)-1]

            s += (math.pow((u_val - uh_val), 2.)*jac*w[j])

    return math.sqrt(s)
