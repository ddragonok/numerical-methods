import numpy as np


def simple_iteration(a, b, n):
    eps = 0.01
    x = np.zeros(n)
    converge = False
    while not converge:
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(a[i][j] * x[j] for j in range(i))
            s2 = sum(a[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / a[i][i]
        converge = np.sqrt(sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_new
    print(x)


def seidel(a, b, n):
    eps = 0.01
    x = np.zeros(n)
    converge = False
    while not converge:
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(a[i][j] * x_new[j] for j in range(i))
            s2 = sum(a[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / a[i][i]
        converge = np.sqrt(sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_new
    print(x)


a_coef = [
    [15, -4, -6, 5],
    [4, -14, -1, 4],
    [7, -7, 27, -8],
    [-3, -3, 2, -14]
    ]
b_coef = [104, 70, 170, 48]
N = 4

simple_iteration(a_coef, b_coef, N)

seidel(a_coef, b_coef, N)
