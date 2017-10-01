import numpy as np
import scipy.sparse as sps


def kab(h, a, b):

    return ((-1)**((a+1)*(b+1))) / h


def fa(h, a, fun, g, e, n):

    if e == 0:
        if a == 0:
            s = 1
        else:
            s = 0
    elif 0 < e < n - 1:
        s = 0
    else:
        s = -kab(h, a, 1) * g

    return (((2.0 - a) * fun[0] + (1.0 + a) * fun[1]) / h) + s


def location_matrix(a, e, n):
    if a == 0:
        return e
    else:
        if e + 1 > n:
            return 0
        else:
            return e + 1


# def map_basis(direction, basis):
#     """
#
#     map_basis
#
#     :param direction: {x_to_c, c_to_x}
#     :type direction: str
#     :return: value of basis in new coordinate
#     """
#
#     if direction == 'x_to_c':
#         xe = 0
#         for a in xrange(2):
#             xe += basis[a](c)*x[a]
#         return xe


def solve_d(n, basis):

    h = 1.0/float(n)
    g = 0
    fun = [1.0, 1.0]
    k = np.empty((n, n), dtype=float)
    f = np.empty((n, 1), dtype=float)

    for e in xrange(n):
        ke = np.empty((2, 2), dtype=float)
        fe = np.empty((2, 1), dtype=float)
        for a in xrange(2):
            for b in xrange(2):
                ke[a][b] = kab(h, a, b)
            fe[a] = fa(h, a, fun, g, e, n)
        for a in xrange(2):
            m = location_matrix(a, e, n)
            if m != 0:
                for b in xrange(2):
                    n = location_matrix(b, e, n)
                    if n != 0:
                        k[m][n] = ke[a][b]
            f[m] = fe[a]

    k = np.asmatrix(k)
    d = k.I * f

    print d

solve_d(4, 4)
