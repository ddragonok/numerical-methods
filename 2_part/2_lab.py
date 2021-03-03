import numpy as np


def tridiagonal_matrix_algorithm(a, b, c, d, n):
    p = np.zeros([n])
    q = np.zeros([n])
    x = np.zeros([n])
    p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]
    for i in range(1, n):
        p[i] = -c[i] / (b[i] + a[i]*p[i - 1])
        q[i] = (d[i] - a[i]*q[i - 1]) / (b[i] + a[i]*p[i - 1])

    x[n - 1] = q[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = p[i]*x[i + 1] + q[i]
    print(x)


a_coef = [0, -3, -2, -4, 4]
b_coef = [12, -18, -16, 18, -9]
c_coef = [-5, -8, -9, -7, 0]
d_coef = [148, 45, -155, 11, 3]
N = 5

tridiagonal_matrix_algorithm(a_coef, b_coef, c_coef, d_coef, N)
