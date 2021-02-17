import math
import numpy as np
import matplotlib.pyplot as plt

# Initial data
a = 2
b = 1
c = -2
N = 10
K = 10
L = math.pi
h = L / N
t = (h*h) / 6
T = t * K

first = (a/(h*h)) + b/h
second = ((1/t) - ((2*a)/(h*h)) - (b/h) + c)
third = a/(h*h)

u = np.zeros((N + 1, K + 1))  # Finite-difference mesh


def create_u():
    for j in range(N + 1):  # Initial condition
        u[j][0] = math.sin(j*h)

    for k in range(K):
        for j in range(1, N):
            u[j][k+1] = (u[j+1][k] * first + u[j][k] * second + u[j-1][k] * third) * t

    for k in range(1, K + 1):  # Boundary conditions
        u[0][k] = (u[1][k] * (1/h) - math.exp((c - a) * k*t) * (math.cos(b * k*t) + math.sin(b * k*t))) / ((1/h) - 1)
        u[N][k] = (u[N-1][k] * (1/h) - math.exp((c - a) * k*t) * (math.cos(b * k*t) + math.sin(b * k*t))) / ((1/h) + 1)

    print(np.around(u, 3))

    x = [h * i for i in range(N + 1)]
    for i in range(K + 1):  # Graph of each time layer
        plt.plot(x, u[:, i])
    plt.show()


create_u()
