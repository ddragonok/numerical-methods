import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
p_x = 0*x
q_x = -2 * (1 + x**2)
f_x = 0*x
a = 0
b = 1
y_diff = 2
y_last = 3 + sym.pi/2
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def gauss(a, b, n):
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


def finite_diff_method(h, a, b, y_diff, y_last, p_x, q_x, f_x):
    x_array = [a, a + h]
    count = int((b - a) / h) + 1
    A = np.zeros((count, count))
    B = np.zeros(count)

    #  First boundary condition
    A[0][0] = -1 / h
    A[0][1] = 1 / h
    B[0] = y_diff

    for i in range(1, count - 1):
        A[i][i - 1] = 1 - (h * p_x.subs(x, x_array[-1])) / 2
        A[i][i] = -2 + h**2 * q_x.subs(x, x_array[-1])
        A[i][i + 1] = 1 + (h * p_x.subs(x, x_array[-1])) / 2
        B[i] = f_x.subs(x, x_array[-1]) * h**2
        x_array.append(x_array[-1] + h)

    #  Second boundary condition
    A[count - 1][count - 1] = 1
    B[count - 1] = y_last

    y_array = gauss(A, B, count)
    return x_array, y_array


diff1 = finite_diff_method(h1, a, b, y_diff, y_last, p_x, q_x, f_x)
diff2 = finite_diff_method(h2, a, b, y_diff, y_last, p_x, q_x, f_x)
draw_graph(diff1, diff2)
