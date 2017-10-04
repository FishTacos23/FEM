import numpy as np
import scipy.sparse as sps


def kab(h, a, b):

    return ((-1)**(a+b)) / h


def fa(h, a, fun, g, hc, e, n, x):

    if e == 1:
        if a == 1:
            s = hc
        else:
            s = 0
    elif 1 < e < n:
        s = 0
    else:
        s = -kab(h, a, 2) * g

    return ((((3 - a) * fun(x)) + (a * fun(x + h))) * (h / 6.0)) + s


def location_matrix(a, e, n):
    if a == 1:
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


def solve_d(n, basis, fun):

    h = 1.0/float(n)
    g = 0
    hc = 0
    k = np.zeros((n, n), dtype=float)
    f = np.zeros((n, 1), dtype=float)
    x = 0

    for e in xrange(1, n+1):
        ke = np.empty((2, 2), dtype=float)
        fe = np.empty((2, 1), dtype=float)
        for a in xrange(1, 3):
            for b in xrange(1, 3):
                ke[a-1][b-1] = kab(h, a, b)
            fe[a-1] = fa(h, a, fun, g, hc, e, n, x)
        for a in xrange(1, 3):
            i = location_matrix(a, e, n)
            if i != 0:
                for b in xrange(1, 3):
                    j = location_matrix(b, e, n)
                    if j != 0:
                        k[i-1][j-1] += ke[a-1][b-1]
                f[i-1] += fe[a-1]
        x += h

    k = np.asmatrix(k)
    d = k.I * f

    print d
