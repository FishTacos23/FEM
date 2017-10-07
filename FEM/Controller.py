import Finite_Element_Method as FEM
import Plotter as Pl
import numpy as np


def func_c(ex):
    return 1.0


def func_x(ex):
    return ex


def func_x2(ex):
    return ex**2


def linear(a, c):

    if a == 1:
        n = 0.5 * (1 - c)

    else:
        n = 0.5 * (1 + c)

    return n


def eq_x(ex, g=0, h=0):
    return (-1. / 6.0) * (ex**3.0 - 1.0) - h*ex + g + h


def eq_x2(ex, g=0, h=0):
    return (-1. / 12.0) * (ex**4.0 - 1.0) - h*ex + g + h


def eq_c(ex, g=0, h=0):
    return (-1. / 2.0) * (ex**2.0 - 1.0) - h*ex + g + h


Model_x2 = FEM.FEM(10, linear, func_c)
uhx2, xhx2 = Model_x2.solve()

x = [np.arange(0., 1., .001)]
u = [eq_x2(x[0])]

Pl.plot_graphs([x, xhx2], [u, uhx2])
