import Finite_Element_Method as FEM
import Plotter as Pl
import numpy as np
import math


def func_p2(ex, h):
    return 10.*(h**3.)


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


def plot_solutions(ps, l, hs, bounds):

    b = 0.005
    mod_e = 10000000.
    row = .1
    pol_i = b * (hs ** 3.) / 12.

    # freq = [float(n)*math.pi*math.sqrt(mod_e/row)/l for n in xrange(1, 1001)]
    freq = [(float(n)-.5)*math.pi*math.sqrt(mod_e/row)/l for n in xrange(1, 1001)]

    n_m_n = np.linspace(1./1000., 1., 1000)

    freq_e_list = []
    n_m_n_list = []

    for j, p in enumerate(ps):
        fp = []
        n_m_n_p = []
        n_adj = 1000-p
        model = FEM.FEM(n_adj, [bern, d_bern, dd_bern], func_p2, l=l, p=p, prop=(pol_i, mod_e, row), h=hs, bc=bounds)
        uh, xh, d, mh, f_list = model.solve()
        # Pl.plot_modes(xh, mh)
        Pl.animation_plot(xh.flatten(), mh[1].flatten())
        for i, f in enumerate(f_list):
            fp.append(f/freq[i])
            n_m_n_p.append(n_m_n[i])
        freq_e_list.append(fp)
        n_m_n_list.append(n_m_n_p)

    Pl.plot_errors(n_m_n_list, freq_e_list)

p_list = [1, 2, 3]
length = 1.
h_list = .005
bc = [[-1, 0]]
# bc = ((0, 0), (-1, 0))

plot_solutions(p_list, length, h_list, bc)
