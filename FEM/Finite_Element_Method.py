import numpy as np
import scipy.sparse as sps


class FEM(object):

    def __init__(self, n, basis, fun, p=1, he=None, h=0, g=0):

        self.n = n
        self.p = p
        self.fun = fun
        self.basis = basis

        if he is None:
            self.he = [1.0/float(n)]*n
        else:
            self.he = he

        self.h = h
        self.g = g

        self.num_basis = n - p - 1
        self.n_int = p + 1

        self.knot_vector = np.array(np.zeros(p+1))
        self.knot_vector = np.insert(self.knot_vector, len(self.knot_vector), np.arange(1, n+1-2*(p+1)))
        self.knot_vector = np.insert(self.knot_vector, len(self.knot_vector), np.zeros(p+1)+(n+1-2*(p+1)))

    def solve(self):

        d = self._solve_d()

        x = []
        u = []

        # x_pos = 0.

        x.append(np.arange(0., 1., 0.000001))
        n_vals = self.get_n(len(x))
        ef = self.xga()

        # for e in xrange(1, self.n + 1):
        #     for j in xrange(1, self.p):
        #         # local B
        #         # local N
        #         # global
        #
        #         for a in xrange(1, self.p+1):
        #             if self._location_matrix(a, e) == 0:
        #                 break
        #             for b in xrange(1, self.p+1):
        #                 if self._location_matrix(b, e) == 0:

            #
            # u.append(ef[e-1]*n_vals[e-1])
            #
            # if e < self.n:
            #     u_e = n_val[0]*d.item(e-1) + n_val[1]*d.item(e)
            # else:
            #     u_e = n_val[0]*d.item(e-1) + n_val[1]*self.g
            #
            # u.append(u_e)
            #
            # x_pos += self.he[e-1]

        print u
        return u, x, d

    def get_n(self, x_len):

        c = np.arange(-1., 1., 2.0 / float(x_len))

        ne = []

        for i in xrange(1, self.n + 1):
            c_mat = np.asmatrix(self.get_c_e(i))
            b_vec = np.empty((self.p+1, 1))
            for j in xrange(1, self.p+1):
                b_vec[j] = self.basis(j, c, self.p)
            ne.append(c_mat*b_vec)

        return ne

    def xga(self):
        x_a = np.empty(self.n-self.p-1)
        for a in xrange(1, self.n-self.p):
            x_a_sum = 0.
            for j in xrange(a+1, self.p+a+1):
                x_a_sum += self.knot_vector[j-1]
            x_a[a-1] = float(x_a_sum)/float(self.p)
        return x_a

    def get_c_e(self, e):

        if self.p == 2:
            if e == 1:
                return np.asarray([[1., 0., 0.], [0., 1., .5], [0., 0., .5]])
            elif e == self.n:
                return np.asarray([[.5, 0., 0.], [.5, 1., 0.], [0., 0., 1.]])
            else:
                return np.asarray([[.5, 0., 0.], [.5, 1., .5], [0., 0., .5]])

        elif self.p == 3:
            if e == 1:
                return np.asarray([[1., 0., 0., 0.], [0., 1., .5, .25],
                                   [0., 0., .5, 7./12.], [0., 0., 0., 1./6.]])
            elif e == 2:
                return np.asarray([[.25, 0., 0., 0.], [7./12., 2./3., 1./3., 1./6.],
                                   [1./6., 1./3., 2./3., 2./3.], [0., 0., 0., 1./6.]])
            elif e == self.n - 1:
                return np.asarray([[1./.6, 0., 0., 0], [2./3., 2./3., 1./3., 1./6.],
                                   [1./6., 1./3., 2./3., 7./12.], [0., 0., 0., .25]])
            elif e == self.n:
                return np.asarray([[1./6., 0., 0., 0.], [7./12., .5, 0., 0.],
                                   [.25, .5, 1., 0.], [0., 0., 0., 1.]])
            else:
                return np.asarray([[1./6., 0., 0., 0.], [2./3., 2./3., 1./3., 1./6.],
                                   [1./6., 1./3., 2./3., 2./3.], [0., 0., 0., 1./6.]])

    def x_to_c(self, ):
        pass

    def _solve_d(self):

        k = np.zeros((self.n, self.n), dtype=float)
        f = np.zeros((self.n, 1), dtype=float)
        x = 0

        for e in xrange(1, self.n + 1):
            for j in xrange(1, self.p):
                # local B
                # local N
                # global

                for a in xrange(1, self.p + 1):
                    i = self._lm(a, e)
                    if i == 0:
                        break
                    for b in xrange(1, self.p + 1):
                        m = self._lm(b, e)
                        if m == 0:
                            k[i - 1][m - 1] += self._kab(a, b, e)
                    f[i - 1] += self._fa(a, e, x)

                    # for e in xrange(1, self.n+1):
        #
        #     for a in xrange(1, 3):
        #         i = self._location_matrix(a, e)
        #         if i != 0:
        #             for b in xrange(1, 3):
        #                 j = self._location_matrix(b, e)
        #                 if j != 0:
        #                     k[i-1][j-1] += self._kab(a, b, e)
        #             f[i-1] += self._fa(a, e, x)
        #
        #     x += self.he[e-1]
        #
        # k = np.asmatrix(k)
        # return k.I * f

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

    def _lm(self, a, e):
        if e + 1 > self.n:
            return 0
        else:
            return e + a - 1
