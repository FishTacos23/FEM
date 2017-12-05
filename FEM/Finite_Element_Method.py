import numpy as np
import math


class FEM(object):

    def __init__(self, n, basis, fun, l=None, p=1, prop=(1., 1., 1., 1., 1., 1., 1.),
                 bc=((-1, 0, 0), (-1, 1, 0), (-1, 2, 0), (-1, 3, 0), (-1, 4, 0), (-1, 5, 0)), h=None):

        self.n = n  # number of elements
        self.p = p  # degree
        self.num_nodes = p + n  # number of nodes
        self.num_basis = n - p - 1  # number of basis functions
        self.n_int = p + 1  # number of quadratures / basis per element
        self.num_bc = len(bc)  # number of boundary conditions

        self.dof = 6  # number of degrees of freedom per node

        self.fun = fun  # forcing function
        self.basis = basis  # basis function
        self.bc = []  # boundary conditions

        for i in xrange(self.num_bc):
            if bc[i][0] >= 0:
                self.bc.append(bc[i])
            else:
                self.bc.append([self.num_nodes+bc[i][0], bc[i][1], bc[i][2]])

        self.i1 = prop[0]  # moment of inertia (I1)
        self.i2 = prop[1]  # moment of inertia (I1)
        self.m_e = prop[2]  # Modulus of Elasticity (E)
        self.v = prop[3]  # poisson's ration (nu)
        self.g = prop[4]  # Shear Modulus of Elasticity (G)
        self.j = prop[5]  # polar moment of inertia (Ip)
        self.a = prop[6]  # Area

        if l is None:
            self.l = float(n)
        else:
            self.l = l  # length of system

        self.knot_vector = self._get_knot_vector()  # knot vector
        self.xga = self._xga()  # global locations of the nodes
        self.id_mat = self._construct_id()  # id array for location matrix
        self.qs, self.ws = self._get_quadratures()  # quadrature points and weighting values
        self.h = h

    def solve(self):

        d_n, ax, n_list = self._solve_d()  # solve for displacements of non-zero nodes

        x = []
        for i in xrange(len(ax)):
            x.append(ax[i][0])

        y = [0.]*len(x)
        z = [0.]*len(x)

        w_theta = []

        d_i = 0
        dn_i = 0
        d = np.zeros((self.num_nodes*self.dof, 1))  # place zero value nodes in d
        for a in xrange(1, len(self.xga)+1):
            for b in xrange(1, self.dof+1):
                n_id = self.id_mat[a-1][b-1]
                if n_id != 0:
                    d[d_i] = d_n[dn_i]
                    dn_i += 1
                d_i += 1

        for e in xrange(1, self.n+1):  # loop over elements
            w_theta_e = np.zeros((len(n_list[e-1]), self.dof))
            for i in xrange(len(n_list[e-1])):  # loop over each xc value
                for ii in xrange(self.dof):
                    for a in xrange(1, self.p + 2):
                        w_theta_e[i][ii] += n_list[e-1][i][a-1]*d[(self._ien(e, a)-1)*self.dof+ii]
            w_theta.append(w_theta_e)

        dx = []
        dy = []
        dz = []

        for e in w_theta:
            for p in e:
                dz.append(p[0]-y[len(dy)]*p[5])
                dy.append(p[1]+x[len(dz)-1]*p[5])
                dx.append(p[2]-x[len(dz)-1]*p[4]+y[len(dy)-1]*p[3])

        return dx, dy, dz, x, y, z

    def _solve_d(self):

        b_d = self._get_big_d()
        k = np.zeros((self.num_nodes*self.dof-self.num_bc, self.num_nodes*self.dof-self.num_bc), dtype=float)
        f = np.zeros((self.num_nodes*self.dof-self.num_bc, 1), dtype=float)  # global F
        x_list = []
        n_list = []

        for e in xrange(1, self.n + 1):  # loop over elements
            n_el_list = []
            ke = np.zeros(((self.p+1)*self.dof, (self.p+1)*self.dof), dtype=float)
            fe = np.zeros(((self.p+1)*self.dof), dtype=float)
            for j in xrange(1, self.n_int + 1):  # loop over quadrature points
                dnx, ddnx, jac, x, ne, dne = self._basis_x(e, self.qs[j - 1])  # get values of global basis, x, and jac
                big_f = self.fun(x)
                x_list.append(x)
                n_el_list.append(ne)
                for a in xrange(1, self.p + 2):  # loop to place element k, f into global K, F
                    ba = self._get_big_b(ne, dnx, a)
                    for b in xrange(1, self.p + 2):  # loop to place element k into global K
                        bb = self._get_big_b(ne, dnx, b)
                        bdb = np.asarray(ba.T * b_d * bb)
                        for i_dof in xrange(1, self.dof+1):
                            for j_dof in xrange(1, self.dof+1):
                                ke[(a-1)*self.dof+i_dof-1][(b-1)*self.dof+j_dof-1] += bdb[i_dof-1][j_dof-1]*jac*self.ws[j-1]
                for a in xrange(1, self.p + 2):
                    for i_dof in xrange(1, self.dof+1):
                            fe[(a-1)*self.dof+i_dof-1] += ne[a-1]*big_f[i_dof-1]*jac*self.ws[j-1]
                for a in xrange(1, self.p+2):
                    for i_dof in xrange(1, self.dof+1):
                        i = self._lm(a, e, i_dof)
                        if i != 0:
                            for b in xrange(1, self.p+2):
                                for j_dof in xrange(1, self.dof+1):
                                    m = self._lm(b, e, j_dof)
                                    if m != 0:
                                        k[i-1][m-1] += ke[(a-1)*self.dof+i_dof-1][(b-1)*self.dof+j_dof-1]
                for a in xrange(1, self.p+2):
                    for i_dof in xrange(1, self.dof+1):
                        i = self._lm(a, e, i_dof)
                        if i != 0:
                            f[i-1] += fe[(a-1)*self.dof+i_dof-1]
            n_list.append(n_el_list)

        k = np.asmatrix(k)
        return np.asarray(k.I * f), x_list, n_list

    def _basis_x(self, e, xc):

        b_s = []  # collection of values for the basis functions at xc
        db_s = []  # collection of values for the derivative of the basis functions at xc
        ddb_s = []  # collection of the values for the double derivative of the basis functions at xc

        for a in xrange(1, self.n_int + 1):  # loop over the local basis
            basis = self._local_b(xc, a)
            b_s.append(basis[0])
            db_s.append(basis[1])
            ddb_s.append(basis[2])

        ce = self._get_c_e(e)
        b_s = np.asarray(b_s)
        db_s = np.asarray(db_s)
        ddb_s = np.asarray(ddb_s)

        ne, dne, ddne = self._local_n(b_s, db_s, ddb_s, ce)  # get local basis

        dnx, ddnx, jac, x = self._global(ne, dne, ddne, self.xga[e - 1:e + self.p + 2])  # convert local to global

        dnx = np.asarray(dnx).flatten()
        ddnx = np.asarray(ddnx).flatten()
        jac = np.asarray(jac)[0][0]
        x = np.asarray(x).flatten()

        return dnx, ddnx, jac, x, ne, dne

    def _get_big_b(self, n, dn, i):
        big_b = np.zeros((self.dof, self.dof), dtype=float)
        big_b[0][2] = dn[i-1]
        big_b[1][0] = dn[i-1]
        big_b[1][4] = -n[i-1]
        big_b[2][1] = dn[i-1]
        big_b[2][3] = n[i-1]
        big_b[3][3] = dn[i-1]
        big_b[4][4] = dn[i-1]
        big_b[5][5] = dn[i-1]
        return np.matrix(big_b)

    def _get_big_d(self):
        big_d = np.zeros((self.dof, self.dof), dtype=float)
        big_d[0][0] = self.m_e*self.a
        big_d[1][1] = self.g*self.a*5./6.
        big_d[2][2] = self.g*self.a*5./6.
        big_d[3][3] = self.m_e*self.i1
        big_d[4][4] = self.m_e*self.i2
        big_d[5][5] = self.g*self.j
        return np.matrix(big_d)

    def _get_knot_vector(self):
        knot_vector = np.array(np.zeros(self.p + 1))
        knot_vector = np.insert(knot_vector, len(knot_vector), np.arange(self.l/self.n, self.l, self.l/self.n))
        knot_vector = np.insert(knot_vector, len(knot_vector), np.zeros(self.p + 1) + self.l)

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
            if self.n == 1:
                return np.asarray([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            elif e == 1:
                return np.asarray([[1., 0., 0.], [0., 1., .5], [0., 0., .5]])
            elif e == self.n:
                return np.asarray([[.5, 0., 0.], [.5, 1., 0.], [0., 0., 1.]])
            else:
                return np.asarray([[.5, 0., 0.], [.5, 1., .5], [0., 0., .5]])

        elif self.p == 3:
            if self.n == 1:
                return np.asarray([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
            elif e == 1:
                return np.asarray([[1., 0., 0., 0.], [0., 1., .5, .25],
                                   [0., 0., .5, 7./12.], [0., 0., 0., 1./6.]])
            elif e == 2:
                return np.asarray([[.25, 0., 0., 0.], [7./12., 2./3., 1./3., 1./6.],
                                   [1./6., 1./3., 2./3., 2./3.], [0., 0., 0., 1./6.]])
            elif e == self.n - 1:
                return np.asarray([[1./6., 0., 0., 0], [2./3., 2./3., 1./3., 1./6.],
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

    def _local_n(self, b, db, ddb, ce):

        ce = np.matrix(ce)
        b = np.matrix(b)
        db = np.matrix(db)
        ddb = np.matrix(ddb)

        n = ce*b.T
        dn = ce*db.T
        ddn = ce*ddb.T

        return np.asarray(n), np.asarray(dn), np.asarray(ddn)

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

        return self.fun(x, self.h)

    def _lm(self, a, e, i_dof):

        global_a = self._ien(e, a)
        return self._id(global_a, i_dof)

    def _ien(self, e, a):
        return e + a - 1

    def _id(self, global_a, i_dof):
        return self.id_mat[global_a - 1][i_dof-1]

    def _construct_id(self):
        id_mat = np.zeros((self.num_nodes, self.dof), dtype=int)
        num_eq = 0
        for i in xrange(self.num_nodes):
            for j in xrange(self.dof):
                skip = False
                for bc in self.bc:
                    if bc[0] == i and bc[1] == j:
                        skip = True
                if not skip:
                    num_eq += 1
                    id_mat[i][j] = num_eq
        return id_mat
