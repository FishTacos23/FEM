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


def bernstein(a, gc, p):
    f1 = (1/(2**float(p)))
    f2 = (math.factorial(p)/(math.factorial(a-1)*math.factorial(p-a+1)))
    f3 = ((1.-gc)**(p-float(a)+1.))
    f4 = ((1.+gc)**(float(a)-1.))

    return f1*float(f2)*f3*f4


def d_bernstein(a, gc, p):
    f1 = (1 / (2 ** float(p)))
    f2 = (math.factorial(p) / (math.factorial(a - 1) * math.factorial(p - a + 1)))
    f3 = ((1. - gc) ** (p - float(a)))
    f4 = ((1. + gc) ** (float(a) - 2.))
    f5 = (2. * float(a) - float(p) * gc - float(p) - 2.)

    return f1*float(f2)*f3*f4*f5


def dd_bernstein(a, gc, p):
    f1 = (1 / (2 ** float(p)))
    f2 = (math.factorial(p) / (math.factorial(a - 1) * math.factorial(p - a + 1)))
    f3 = ((1. - gc) ** (p - float(a) - 1.))
    f4 = ((1. + gc) ** (float(a) - 3.))
    f5 = (4.*(float(a)**2.)-4.*float(a)*(float(p)*(gc+1.)-gc+2.) +
          (float(p)**2.)*((gc + 1.)**2.)+float(p)*(-gc-1.)*(gc-3.)-4.*(gc-1.))

    return f1 * float(f2) * f3 * f4 * f5


def eq_x2(ex, g=0, h=0):
    return (-1. / 12.0) * (ex**4.0 - 1.0) - h*ex + g + h


def plot_errors(ns, funcs, eqs, ps, l):

    hs = []
    dhs = []
    error = []

    for n in ns:
        hs.append(l/float(n))

    for func in funcs:
        for p in ps:
            e_f = []
            dh_f = []
            for n in ns:
                model = FEM.FEM(n, [bernstein, d_bernstein, dd_bernstein], func, l=l, p=p)
                uh, xh, d, xga = model.solve()
                e_f.append(Er.calc_error(None, d, linear, eqs[0], model, p=p))
                dh_f.append(d)
            dhs.append(dh_f)
            error.append(e_f)

    Pl.plt_error(hs, error, 'Error')
    Pl.plt_error(ns, error, 'Error')


def plot_solutions(ns, funcs, eqs, ps, l):

    x = [np.arange(0., 1., .000001)]

    for i, func in enumerate(funcs):
        for p in ps:
            for n in ns:
                model = FEM.FEM(n, [bernstein, d_bernstein, dd_bernstein], func, l=l, p=p)
                uh, xh, d, xga = model.solve()

                Pl.plot_graphs([x, xh], [[eqs[i](x[0])], uh], 'n=' + str(n) + ' f=x2', p=[d, xga])


n_list = [1, 10, 100, 1000]
func_list = [func_x2]
eq_list = [eq_x2]
p_list = [2, 3]
length = 1.

plot_errors(n_list, func_list, eq_list, p_list, length)
# plot_solutions(n_list, func_list, eq_list, p_list, length)
