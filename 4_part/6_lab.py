import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

# Initial data
a = 0
b = 1
x = sym.symbols('x')
t = sym.symbols('t')
K = x * t
f = (5/6) * x
n = 50
h = (b - a) / n

u_app = np.zeros(n + 1)


def draw_u(x, u):
    plt.grid(True)
    plt.plot(x, u)
    plt.show()


def create_u(u, K, f):
    x_array = [(h*i + a) for i in range(n + 1)]
    wj = 1
    A = np.zeros((n + 1, n + 1))
    B = np.zeros(n + 1)

    for i in range(n + 1):
        for j in range(n + 1):
            A[i][j] = -h * wj * K.subs([(x, x_array[i]), (t, x_array[j])])
        A[i][i] = A[i][i] + 1
    for j in range(n + 1):
        B[j] = f.subs(x, x_array[j])
    A_inverse = np.linalg.inv(A)
    u = B.dot(A_inverse)

    print(np.around(u, 3))
    draw_u(x_array, u)


create_u(u_app, K, f)