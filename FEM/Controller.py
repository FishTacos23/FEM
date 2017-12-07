import Finite_Element_Method as FEM
import Plotter as Pl
import numpy as np
import Error as Er
import math


# TODO questions: Can I hard code 6 dof? Can I look for


def func_x2(ex):
    return ex**2.


def func_p2(ex, h):
    return 10.*(h**3.)


def linear(a, c, p):

    if a == 1:
        n = 0.5 * (1 - c)

    else:
        n = 0.5 * (1 + c)

    return n


def uniform_axial_loading(ex):
    f = np.asarray([0., 0., 100000., 0., 0., 0.])
    return f


def non_uniform_axial_loading(ex):
    f = np.asarray([0., 0., 100000.*ex, 0., 0., 0.])
    return f


def transverse_loading(ex):
    return np.asarray([0., 10000., 0., 0., 0., 0.])


def uniform_axial_deformation(ex):
    return 100000.


def non_uniform_axial_deformation(ex):
    return 100000.*ex


def transverse_deformation(ex, e_, i_):
    return 10000.*(ex**4.)/(8.*e_*i_)


def transverse_rot(ex, e_, i_):
    return 10000.*(ex**3.)/(6.*e_*i_)


def bern(a, gc, p):
    f1 = (1/(2**float(p)))
    f2 = (math.factorial(p)/(math.factorial(a-1)*math.factorial(p-a+1)))
    f3 = ((1.-gc)**(p-float(a)+1.))
    f4 = ((1.+gc)**(float(a)-1.))

    return f1*float(f2)*f3*f4


def d_bern(a, gc, p):
    f1 = (1 / (2 ** float(p)))
    f2 = (math.factorial(p) / (math.factorial(a - 1) * math.factorial(p - a + 1)))
    f3 = ((1. - gc) ** (p - float(a)))
    f4 = ((1. + gc) ** (float(a) - 2.))
    f5 = (2. * float(a) - float(p) * gc - float(p) - 2.)

    return f1*float(f2)*f3*f4*f5


def dd_bern(a, gc, p):
    f1 = (1 / (2 ** float(p)))
    f2 = (math.factorial(p) / (math.factorial(a - 1) * math.factorial(p - a + 1)))
    f3 = ((1. - gc) ** (p - float(a) - 1.))
    f4 = ((1. + gc) ** (float(a) - 3.))
    f5 = (4.*(float(a)**2.)-4.*float(a)*(float(p)*(gc+1.)-gc+2.) +
          (float(p)**2.)*((gc + 1.)**2.)+float(p)*(-gc-1.)*(gc-3.)-4.*(gc-1.))

    return f1 * float(f2) * f3 * f4 * f5


def eq_x2(ex, g=0, h=0):
    return (-1. / 12.0) * (ex**4.0 - 1.0) - h*ex + g + h


def beam_theory(ex, h):
    w = func_p2(ex, h)
    b = 0.005
    mod_e = 1000000.
    pol_i = b * (h ** 3.) / 12.
    l = 1.

    return w*(ex**4 - 4*(l**3)*ex + 3*l**4)/(24*mod_e*pol_i)


def axial_beam_theory(ex, e, a):
    n = uniform_axial_deformation(ex)
    return n*ex/(e*a)


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


def plot_solutions(ns, funcs, eqs, ps, l, hs, bc):

    x = [np.arange(0., l, .0001)]
    max_deflection = []
    b = 0.005
    mod_e = 10000000.
    v = .33
    g = 3770000.
    r = .1
    a = math.pi*(r**2.)
    i1 = math.pi*(r**4.)/4.
    i2 = i1
    ip = i1*2.
    exact = axial_beam_theory(l, mod_e, a)
    exact_def = [transverse_deformation(xi, mod_e, i1) for xi in x]
    exact_rot = [transverse_rot(xi, mod_e, i1) for xi in x]
    defs = []
    rots = []
    labels = []
    x_s = []

    x_s.append(x)
    defs.append(exact_def)
    rots.append(exact_rot)
    labels.append('EXACT SOLUTION')

    for h in hs:
        for i, func in enumerate(funcs):
            for j, p in enumerate(ps):
                max_n_deflect = []
                for n in ns:

                    model = FEM.FEM(n, [bern, d_bern, dd_bern], func, l=l, p=p,
                                    prop=(i1, i2, mod_e, v, g, ip, a), h=h, bc=bc)
                    dx, dy, dz, x, y, z, xt, yt, zt = model.solve()
                    max_n_deflect.append(max(dy))
                    # Pl.plot_graphs([x, xh], [[eqs[i](x[0], h)], uh], 'n=' + str(n) + ' f=x2', p=[d, xga])
                    # Pl.animated_plot(x, dx, y, dy)

                    # p_x = []
                    # for q in xrange(len(x)):
                    #     p_x.append(x[q]+dx[q])
                    # defs.append(p_x)

                    # defs.append(dy)
                    # rots.append(zt)
                    # x_s.append(x)
                    # labels.append('P='+str(p)+' N='+str(n))
                max_deflection.append(max_n_deflect)

    Pl.plot_deflections(ns, max_deflection, exact_def[0][-1], 'Deflection at Tip of Cantilever Beam', ['N=1', 'N=10', 'Exact'])
    # Pl.new_plot(x_s, defs, rots, labels)

# func_list = [uniform_axial_loading]
func_list = [transverse_loading]
eq_list = [beam_theory]
# p_list = [1]
# p_list = [2]
p_list = [1, 2, 3]
length = 1.
# n_list = [1, 10]
n_list = [10, 100, 1000]
# n_list = [10]
h_list = [.005]
bc = ((0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0), (0, 5, 0),)
# bc = ((0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 5, 0), (-1, 0, 0), (-1, 1, 0))

plot_solutions(n_list, func_list, eq_list, p_list, length, h_list, bc)
