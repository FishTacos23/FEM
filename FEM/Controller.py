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


def plot_errors(ns, funcs, eqs):

    dhs = []

    for func in funcs:
        dh_f = []
        for n in ns:
            model = FEM.FEM(n, linear, func)
            uh, xh, dh = model.solve()
            dh_f.append(dh)
        dhs.append(dh_f)

    e = []
    for i in xrange(len(funcs)):
        e_f = []
        for j in xrange(len(ns)):
            he = [1.0 / float(ns[j])] * ns[j]
            e_f.append(Er.calc_error(he, dhs[i][j], linear, eqs[i]))
        e.append(e_f)

    hs = []
    for n in ns:
        hs.append(1. / n)

    for e_func in e:
        for e_n in e_func:
            print e_n

    Pl.plt_error(hs, e, 'Error')


def plot_solutions(ns, funcs, eqs):

    x = [np.arange(0., 1., .000001)]

    for i, func in enumerate(funcs):
        for n in ns:

            model = FEM.FEM(n, linear, func)
            uh, xh, dh = model.solve()

            Pl.plot_graphs([x, xh], [[eqs[i](x[0])], uh], 'n=' + str(n) + ' f=x2')



n_list = [10, 100, 1000, 10000]
func_list = [func_c, func_x, func_x2]
eq_list = [eq_c, eq_x, eq_x2]

plot_errors(n_list, func_list, eq_list)
