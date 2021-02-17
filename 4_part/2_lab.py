import math
import numpy as np
import matplotlib.pyplot as plt

# Initial data
N = 10
K = 10
L = math.pi
h = L / N
t = h / 2
T = t * K

first = (1/(h*h)) + 1/h
second = (2/(t*t)) + (3/t) - (2/(h*h)) - (1/h) - 1
div = (1/(t*t)) + (3/t)

u = np.zeros((N + 1, K + 1))  # Finite-difference mesh


def create_u():
    for j in range(1, N):  # Initial conditions
        u[j][0] = math.sin(j*h)
        u[j][1] = u[j][0] - t*math.sin(j*h)

    for k in range(1, K):
        for j in range(1, N):
            u[j][k+1] = (u[j+1][k]*first + u[j][k]*second + (u[j-1][k]/(h*h)) - (u[j][k-1]/(t*t)) - math.cos(j*h)*math.exp(-(k*t))) / div

    for k in range(K + 1):  # Boundary conditions
        u[0][k] = u[1][k] - h*math.exp(-(k*t))
        u[N][k] = u[N-1][k] - h*math.exp(-(k*t))

    print(np.around(u, 3))

    x = [h * i for i in range(N + 1)]
    for i in range(K + 1):  # Graph of each time layer
        plt.plot(x, u[:, i])
    plt.show()


create_u()
