import math
import numpy as np
import matplotlib.pyplot as plt

# Initial data
a = 3
b = 1
c = 1
d = -1
N = 10
K = 10
L = math.pi
h = L / N
t = h / 2
T = t * K

first = (b/(h*h)) + (c/h)
second = (2/(t*t)) + (a/t) - ((2*b)/(h*h)) - (c/h) + d
div = (1/(t*t)) + (a/t)

u_mesh = np.zeros((N + 1, K + 1))  # Finite-difference mesh


def draw_u(u):
    x = [h * i for i in range(N + 1)]
    for i in range(K + 1):  # Graph of each time layer
        plt.plot(x, u[:, i])
    plt.show()


def create_u(u):
    for j in range(1, N):  # Initial conditions
        u[j][0] = math.sin(j*h)
        u[j][1] = u[j][0] - t*math.sin(j*h)

    for k in range(1, K):
        for j in range(1, N):
            u[j][k+1] = (u[j+1][k]*first + u[j][k]*second + ((u[j-1][k]*b)/(h*h)) - (u[j][k-1]/(t*t)) - math.cos(j*h)*math.exp(-(k*t))) / div

    for k in range(K + 1):  # Boundary conditions
        u[0][k] = u[1][k] - h*math.exp(-(k*t))
        u[N][k] = u[N-1][k] - h*math.exp(-(k*t))

    print(np.around(u, 3))
    draw_u(u)


create_u(u_mesh)
