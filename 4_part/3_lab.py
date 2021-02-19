import math
import numpy as np
import matplotlib.pyplot as plt

# Initial data
eps = 0.001
N_x = 10
N_y = 10
l_x = math.pi/2
l_y = math.pi/2
h_x = l_x / N_x
h_y = l_y / N_y

first = (1/(h_x*h_x)) + (2/h_x)
second = (1/(h_y*h_y)) + (2/h_y)
dev = (2/(h_x*h_x)) + (2/(h_y*h_y)) + (2/h_x) + (2/h_y) - 4

u = np.zeros((N_x + 1, N_y + 1))  # Finite-difference mesh


def create_u():
    for i in range(N_x + 1):  # Boundary conditions for x
        u[i][0] = math.exp(-(i*h_x)) * math.cos(i*h_x)
        u[i][N_y] = 0
    for j in range(1, N_y):  # Boundary conditions for y
        u[0][j] = math.exp(-(j*h_y)) * math.cos(j*h_y)
        u[N_x][j] = 0

    # Simple-iteration method
    u_prev = u.copy()
    for i in range(1, N_x):
        for j in range(1, N_y):
            u_prev[i][j] = 1.0  # Initial approximation

    while True:
        for i in range(1, N_x):
            for j in range(1, N_y):
                u[i][j] = (u_prev[i+1][j]*first + (u_prev[i-1][j]/(h_x*h_x)) + u_prev[i][j+1]*second + (u_prev[i][j-1]/(h_y*h_y))) / dev

        max_u = 0
        for i in range(1, N_x):
            for j in range(1, N_y):
                if max_u < abs(u[i][j] - u_prev[i][j]):
                    max_u = abs(u[i][j] - u_prev[i][j])
        if max_u <= eps:
            break
        else:
            u_prev = u.copy()

    print(np.around(u, 3))

    x = [h_x * i for i in range(N_x + 1)]
    for i in range(N_y + 1):  # Graph of each y-layer
        plt.plot(x, u[:, i])
    plt.show()


create_u()