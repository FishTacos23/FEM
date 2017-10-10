import numpy as np
import scipy.sparse as sps


class FEM(object):

    def __init__(self, n, basis, fun, he=None, h=0, g=0):

        self.n = n
        self.fun = fun
        self.basis = basis

        if he is None:
            self.he = [1.0/float(n)]*n
        else:
            self.he = he

        self.h = h
        self.g = g

    def solve(self):

        d = self._solve_d()

        x = []
        u = []

        x_pos = 0.

        for e in xrange(1, self.n + 1):

            x.append(np.arange(x_pos, x_pos + self.he[e-1], 0.000001))

            n_val = self.get_n(x[-1])

            if e < self.n:
                u_e = n_val[0]*d.item(e-1) + n_val[1]*d.item(e)
            else:
                u_e = n_val[0]*d.item(e-1) + n_val[1]*self.g

            u.append(u_e)

            x_pos += self.he[e-1]

        return u, x, d

    def get_n(self, x):

        c = np.arange(-1., 1., 2.0 / float(len(x)))

        n2 = self.basis(2, c)
        n1 = self.basis(1, c)

        return [n1, n2]

    def x_to_c(self, ):
        pass

    def _solve_d(self):

        k = np.zeros((self.n, self.n), dtype=float)
        f = np.zeros((self.n, 1), dtype=float)
        x = 0

        for e in xrange(1, self.n+1):

            for a in xrange(1, 3):
                i = self._location_matrix(a, e)
                if i != 0:
                    for b in xrange(1, 3):
                        j = self._location_matrix(b, e)
                        if j != 0:
                            k[i-1][j-1] += self._kab(a, b, e)
                    f[i-1] += self._fa(a, e, x)

            x += self.he[e-1]

        k = np.asmatrix(k)
        return k.I * f

    def _kab(self, a, b, e):

        return ((-1)**(a+b)) / self.he[e-1]

    def _fa(self, a, e, x):

        if e == 1:
            if a == 1:
                s = self.h
            else:
                s = 0
        elif 1 < e < self.n:
            s = 0
        else:
            s = -self._kab(a, 2, e) * self.g

        return ((((3 - a) * self.fun(x)) + (a * self.fun(x + self.he[e-2]))) * (self.he[e-1] / 6.0)) + s

    def _location_matrix(self, a, e):
        if a == 1:
            return e
        else:
            if e + 1 > self.n:
                return 0
            else:
                return e + 1
