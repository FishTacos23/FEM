import Finite_Element_Method as FEM
import Plotter as Pl
import numpy as np
import Error as Er


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


ns = [10, 100, 1000, 10000]
funcs = [func_c, func_x, func_x2]
eqs = [eq_c, eq_x, eq_x2]
uhs = []
xhs = []
dhs = []

for func in funcs:
    uh_f = []
    xh_f = []
    dh_f = []
    for n_ in ns:
        Model = FEM.FEM(n_, linear, func)
        uh, xh, dh = Model.solve()
        uh_f.append(uh)
        xh_f.append(xh)
        dh_f.append(dh)
    uhs.append(uh_f)
    xhs.append(xh_f)
    dhs.append(dh_f)

x = [np.arange(0., 1., .000001)]
ucs = [[eq_c(x[0])], [eq_x(x[0])], [eq_x2(x[0])]]

# Pl.plot_graphs([x, xh], [u, uh], 'u='+str(n)+' f=x2')

e = []
for i in xrange(len(funcs)):
    e_f = []
    for j in xrange(len(ns)):
        e_f.append(Er.calc_error(uhs[i][j], xhs[i][j], dhs[i][j], linear, eqs[i]))
    e.append(e_f)

hs = []
for n_1 in ns:
    hs.append(1./n_1)

labels = ['c', 'x', 'x2']

Pl.plt_error(hs, e, labels, 'Error')
print e
