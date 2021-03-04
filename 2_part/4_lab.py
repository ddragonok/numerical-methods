import sympy as sym

a = 0
b = 2
x = sym.Symbol('x')
f = sym.exp(x) - 2*x - 2
fi = sym.log(2*x + 2)


def half_division(f, a, b):
    eps = 0.001
    while abs(b - a) > eps:
        c = (a + b) / 2
        if f.subs(x, a) * f.subs(x, c) < 0:
            b = c
        else:
            a = c
    root = (a + b) / 2
    print(root)


def simple_iteration(fi, a, b):
    eps = 0.001
    x0 = (a + b) / 2
    x_next = fi.subs(x, x0)
    while abs(x_next - x0) > eps:
        x0 = x_next
        x_next = fi.subs(x, x0)
    print(x_next)


def newton(f):
    eps = 0.001
    x0 = 1.9
    f_diff = sym.diff(f)
    x_next = x0 - (f.subs(x, x0) / f_diff.subs(x, x0))
    while abs(x_next - x0) > eps:
        x0 = x_next
        x_next = x0 - (f.subs(x, x0) / f_diff.subs(x, x0))
    print(x_next)


def secant(f):
    eps = 0.001
    x0 = 0.5
    x1 = 1.9
    x_next = x1 - (f.subs(x, x1) / (f.subs(x, x1) - f.subs(x, x0))) * (x1 - x0)
    while abs(x_next - x0) > eps:
        x0 = x1
        x1 = x_next
        x_next = x1 - (f.subs(x, x1) / (f.subs(x, x1) - f.subs(x, x0))) * (x1 - x0)
    print(x_next)


half_division(f, a, b)
simple_iteration(fi, a, b)
newton(f)
secant(f)
