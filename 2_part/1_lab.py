import numpy as np


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

    print(x)


N = 4
a_matrix = [
    [-5, -6, 4, -2],
    [0, 3, -4, -6],
    [2, 4, -4, 2],
    [1, -8, 2, 8]
    ]
b_matrix = [64, -55, -48, 68]

gauss(a_matrix, b_matrix, N)
