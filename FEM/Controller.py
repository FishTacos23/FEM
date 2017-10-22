import Finite_Element_Method as FEM
import Plotter as Pl
import numpy as np
import Error as Er
import math


def func_x2(ex):
    return ex**2


def linear(a, c, p):

    if a == 1:
        n = 0.5 * (1 - c)

    else:
        n = 0.5 * (1 + c)

    return n


def b_spline(p):

    if p == 1:
        pass
    elif p == 2:
        pass
    elif p == 3:
        pass
    else:
        pass

    return 0


def bernstein(a, gc, p):
    f1 = (1/(2**float(p)))
    f2 = (math.factorial(p)/(math.factorial(a-1)*math.factorial(p-a+1)))
    f3 = ((1.-gc)**(p-float(a)+1.))
    f4 = ((1.+gc)**(float(a)-1.))

    return f1*float(f2)*f3*f4


def eq_x2(ex, g=0, h=0):
    return (-1. / 12.0) * (ex**4.0 - 1.0) - h*ex + g + h


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


def plot_solutions(ns, funcs, eqs, ps):

    x = [np.arange(0., 1., .000001)]

    for i, func in enumerate(funcs):
        for n in ns:
            for p in ps:
                model = FEM.FEM(n, bernstein, func, p=p)
                uh, xh, dh = model.solve()

                Pl.plot_graphs([x, xh], [[eqs[i](x[0])], uh], 'n=' + str(n) + ' f=x2')


n_list = [10]
func_list = [func_x2]
eq_list = [eq_x2]
p_list = [1]

plot_solutions(n_list, func_list, eq_list, p_list)
