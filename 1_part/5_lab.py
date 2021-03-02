import numpy as np
import matplotlib.pyplot as plt

inp_x = [-1, 0, 1, 2, 3, 4]
inp_y = [-0.5, 0, 0.5, 0.86603, 1, 0.86603]
N = len(inp_x)


def draw_graph(n, a_coef):
    count = 100
    x = []
    y = []
    f = 0
    h = (inp_x[N - 1] - inp_x[0]) / count
    for i in range(count):
        x.append(inp_x[0] + h*i)
        for j in range(n + 1):
            f += a_coef[j]*pow(x[-1], j)
        y.append(f)
    plt.plot(x, y)
    plt.show()


def gauss(n, a, b):
    eps = 0.0001
    a_coef = np.zeros(n + 1)
    k = 0
    while k <= n:
        for i in range(k, n + 1):
            temp = a[i][k]
            if abs(temp) < eps:
                continue
            for j in range(n + 1):
                a[i][j] = a[i][j] / temp
            b[i] = b[i] / temp
            if i == k:
                continue
            for j in range(n + 1):
                a[i][j] = a[i][j] - a[k][j]
            b[i] = b[i] - b[k]
        k += 1

    for k in range(n, -1, -1):
        a_coef[k] = b[k]
        for i in range(k):
            b[i] = b[i] - a[i][k]*a_coef[k]

    print(a_coef)

    draw_graph(n, a_coef)


def create_matrix(x, y, n):
    a = np.zeros((n + 1, n + 1))
    b = np.zeros(n + 1)
    summa1 = 0
    summa2 = 0
    for i in range(n + 1):
        for k in range(n + 1):
            for j in range(N):
                summa1 += pow(x[j], k + i)
                summa2 += y[j] * pow(x[j], k)
            a[i][k] = summa1
            if i == 0:
                b[k] = summa2
    print('Матрица коэффициентов для n =', n)
    print(a, '\n', b)
    gauss(n, a, b)


create_matrix(inp_x, inp_y, 1)
create_matrix(inp_x, inp_y, 2)
create_matrix(inp_x, inp_y, 3)
