import sympy as sym
import numpy as np

x1 = sym.Symbol('x1')
x2 = sym.Symbol('x2')
f1 = x1 - sym.cos(x2) - 3
f2 = x2 - sym.sin(x1) - 3
fi1 = sym.cos(x2) + 3
fi2 = sym.sin(x1) + 3


def simple_iteration(fi1, fi2):
    eps = 0.001
    x = [0.01, 0.01]
    x_next = [0, 0]
    x_next[0] = fi1.subs(x2, x[1])
    x_next[1] = fi2.subs(x1, x[0])
    while max(abs(x_next[i] - x[i]) for i in range(len(x))) > eps:
        x[0] = x_next[0]
        x[1] = x_next[1]
        x_next[0] = fi1.subs(x2, x[1])
        x_next[1] = fi2.subs(x1, x[0])
    print(x_next)


def gauss(x0, f1, f2):
    diff_11 = sym.diff(f1, x1)
    diff_12 = sym.diff(f1, x2)
    diff_21 = sym.diff(f2, x1)
    diff_22 = sym.diff(f2, x2)
    a = [[diff_11.subs([(x1, x0[0]), (x2, x0[1])]), diff_12.subs([(x1, x0[0]), (x2, x0[1])])],
         [diff_21.subs([(x1, x0[0]), (x2, x0[1])]), diff_22.subs([(x1, x0[0]), (x2, x0[1])])]]
    b = [-f1.subs([(x1, x0[0]), (x2, x0[1])]), -f2.subs([(x1, x0[0]), (x2, x0[1])])]
    n = len(b)
    eps = 0.0001
    x = np.zeros(n)
    k = 0

    while k <= n:
        for i in range(k, n):
            temp = a[i][k]
            if abs(temp) < eps:
                continue
            for j in range(n):
                a[i][j] = a[i][j] / temp
            b[i] = b[i] / temp
            if i == k:
                continue
            for j in range(n):
                a[i][j] = a[i][j] - a[k][j]
            b[i] = b[i] - b[k]
        k += 1

    for k in range(n - 1, -1, -1):
        x[k] = b[k]
        for i in range(k):
            b[i] = b[i] - a[i][k]*x[k]

    return x


def newton(f1, f2):
    eps = 0.001
    x = [1.0, 1.0]
    x_next = [0, 0]
    delta = gauss(x, f1, f2)
    x_next[0] = x[0] + delta[0]
    x_next[1] = x[1] + delta[1]
    while max(abs(x_next[i] - x[i]) for i in range(len(x))) > eps:
        x[0] = x_next[0]
        x[1] = x_next[1]
        delta = gauss(x, f1, f2)
        x_next[0] = x[0] + delta[0]
        x_next[1] = x[1] + delta[1]
    print(x_next)


simple_iteration(fi1, fi2)
newton(f1, f2)
