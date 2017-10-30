import numpy as np
import math
import matplotlib.pyplot as plt


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

        self.geo_i = 1
        self.mat_e = 1

        self.num_basis = n - p - 1
        self.n_int = p + 1

        self.knot_vector = self._get_knot_vector()

    def solve(self):

        d, n = self._solve_d()

        u = []
        xc = np.arange(0., 1., .00001)

        for a in xrange(1, len(d)):
            u_h

        return u, d, x

    def get_n(self, x_len):

        c = np.arange(-1., 1., 2.0 / float(x_len))

        ne = []

        for i in xrange(1, self.n + 1):
            c_mat = np.asmatrix(self._get_c_e(i))
            b_vec = np.empty((self.p+1, 1))
            for j in xrange(1, self.p+1):
                b_vec[j] = self.basis(j, c, self.p)
            ne.append(c_mat*b_vec)

        return ne

    def _solve_d(self):

        xa = self._xga()

        k = np.zeros((len(xa) - 2, len(xa) - 2), dtype=float)
        f = np.zeros((len(xa) - 2, 1), dtype=float)

        q_s, w_s = self._get_quadratures()

        xc_values = np.linspace(-1., 1., 10000)
        b_local_curves = []
        for a in xrange(1, self.n_int + 1):
            b_local_curves.append(self._local_b_curves(xc_values, a))

        for e in xrange(1, self.n+1):
            ce = self._get_c_e(e)


        n_curves = []

        for e in xrange(1, self.n + 1):

            for j in xrange(1, self.n_int + 1):

                b_s = []
                db_s = []
                ddb_s = []

                for a in xrange(1, self.n_int + 1):
                    basis = self._local_b(q_s[j-1], a)
                    b_s.append(basis[0])
                    db_s.append(basis[1])
                    ddb_s.append(basis[2])

                ce = self._get_c_e(e)
                b_s = np.asarray(b_s)
                db_s = np.asarray(db_s)
                ddb_s = np.asarray(ddb_s)

                ne, dne, ddne, n = self._local_n(b_s, db_s, ddb_s, ce)

                dnx, ddnx, jac, x = self._global(ne, dne, ddne, xa[e-1:e+self.p+2])

                dnx = np.asarray(dnx)
                ddnx = np.asarray(ddnx)
                jac = np.asarray(jac)
                x = np.asarray(x)

                for a in xrange(1, self.p + 1):
                    i = self._lm(a, e)
                    if i == 0:
                        break
                    for b in xrange(1, self.p + 1):
                        m = self._lm(b, e)
                        if m == 0:
                            break
                        k[i-1][m-1] = ddnx[a-1][0]*self.mat_e*self.geo_i*ddnx[b-1][0]*jac[0][0]*w_s[j-1]
                    f[i - 1] += ne*self._fa(x[0])*jac[0][0]*w_s[j-1]

        k = np.asmatrix(k)
        return np.asarray(k.I * f), n_curves

    def _get_knot_vector(self):
        knot_vector = np.array(np.zeros(self.p + 1))
        knot_vector = np.insert(knot_vector, len(knot_vector), np.arange(1, self.n))
        knot_vector = np.insert(knot_vector, len(knot_vector), np.zeros(self.p + 1) + self.n)

        return knot_vector

    def _get_quadratures(self):

        if self.n_int == 1:
            q_points = [0.]
            weights = [2]
        elif self.n_int == 2:
            q_points = [-1./math.sqrt(3.), 1./math.sqrt(3.)]
            weights = [1, 1]
        elif self.n_int == 3:
            q_points = [-math.sqrt(3./5.), 0., math.sqrt(3./5.)]
            weights = [5./9., 8./9., 5./9.]
        elif self.n_int == 4:
            q_points = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]
            weights = [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]
        else:
            q_points = None
            weights = None

        return q_points, weights

    def _xga(self):
        x_a = np.empty(self.n+self.p)
        for a in xrange(1, self.n+self.p+1):
            x_a_sum = 0.
            for j in xrange(a+1, self.p+a+1):
                x_a_sum += self.knot_vector[j-1]
            x_a[a-1] = float(x_a_sum)/float(self.p)
        return x_a

    def _get_c_e(self, e):

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

    def _local_b(self, gc, a):

        b = self.basis[0](a, gc, self.p)
        db = self.basis[1](a, gc, self.p)
        ddb = self.basis[2](a, gc, self.p)

        return b, db, ddb

    def _local_b_curves(self, xs, a):
        b_curves = np.empty((len(xs)))
        for i, x in enumerate(xs):
            b_curves[i] = self.basis[0](a, x, self.p)
        return b_curves

    def _local_n(self, b, db, ddb, ce):

        ce = np.matrix(ce)
        b = np.matrix(b)
        db = np.matrix(db)
        ddb = np.matrix(ddb)

        n = ce*b.T
        dn = ce*db.T
        ddn = ce*ddb.T

        return np.asarray(n), np.asarray(dn), np.asarray(ddn), n

    def _global(self, n, dn, ddn, xae):

        xe_c = 0.
        dxe_c = 0.
        ddxe_c = 0.

        for i in xrange(1, self.p+2):
            xe_c += xae[i-1]*n[i-1]
            dxe_c += xae[i-1]*dn[i-1]
            ddxe_c += xae[i-1]*ddn[i-1]

        dxe_c_s = dxe_c**2
        dxe_c_s = np.matrix(dxe_c_s)

        dn = np.matrix(dn)
        dxe_c = np.matrix(dxe_c)

        ddn = np.matrix(ddn)
        ddxe_c = np.matrix(ddxe_c)

        dnx = dn*dxe_c.I
        ddnx = (ddn - dnx*ddxe_c)*dxe_c_s.I
        j = dxe_c

        return dnx, ddnx, j, xe_c

    def _fa(self,  x):

        return self.fun(x)

    def _lm(self, a, e):

        if e + 1 > self.n:
            return 0
        else:
            return e + a - 1
