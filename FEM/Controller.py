import Finite_Element_Method as FEM


def func_c(x):
    return 1.0


def func_x(x):
    return x


def func_x2(x):
    return x**2


FEM.solve_d(10000, None, func_x2)
