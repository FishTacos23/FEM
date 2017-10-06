import Finite_Element_Method as FEM
import Plotter as Pl


def func_c(x):
    return 1.0


def func_x(x):
    return x


def func_x2(x):
    return x**2


def linear(a, c):

    if a == 1:
        n = 0.5 * (1 - c)

    else:
        n = 0.5 * (1 + c)

    return n

Model = FEM.FEM(4, None, func_x)
u = Model.solve()
Pl.print_value(u)
