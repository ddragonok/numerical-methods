import math
import numpy as np
import matplotlib.pyplot as plt

# Initial data
a = -2
b = -2
c = -4
eps = 0.001
N_x = 10
N_y = 10
l_x = math.pi/2
l_y = math.pi/2
h_x = l_x / N_x
h_y = l_y / N_y

first = (1/(h_x*h_x)) - (a/h_x)
second = (1/(h_y*h_y)) - (b/h_y)
dev = (2/(h_x*h_x)) + (2/(h_y*h_y)) - (a/h_x) - (b/h_y) + c

u_mesh = np.zeros((N_x + 1, N_y + 1))  # Finite-difference mesh


def draw_u(u):
    x = [h_x * i for i in range(N_x + 1)] * (N_y + 1)
    y = np.zeros((N_x + 1) * (N_y + 1))
    i = 0
    j = 0
    while i < (N_x + 1) * (N_y + 1):
        y[i] = h_y * j
        i += 1
        if i % (N_x + 1) == 0:
            j += 1

    z = np.zeros((N_x + 1) * (N_y + 1))
    k = 0
    n = 0
    for i in range((N_x + 1) * (N_y + 1)):
        z[i] = u[n][k]
        n += 1
        if n % (N_x + 1) == 0:
            n = 0
            k += 1

    x2 = np.zeros((N_x + 1) * (N_y + 1))
    i = 0
    j = 0
    while i < (N_x + 1) * (N_y + 1):
        x2[i] = h_x * j
        i += 1
        if i % (N_y + 1) == 0:
            j += 1

    y2 = [h_y * i for i in range(N_y + 1)] * (N_x + 1)

    z2 = np.zeros((N_x + 1) * (N_y + 1))
    i = 0
    j = 0
    r = 0
    while r < (N_x + 1) * (N_y + 1):
        z2[r] = z[i * (N_x + 1) + j]
        r += 1
        i += 1
        if r % (N_y + 1) == 0:
            i = 0
            j += 1

    ax = plt.axes(projection="3d")
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u')
    for i in range(N_y + 1):
        ax.plot3D(x[(N_x + 1) * i:(N_x + 1) * (i + 1)], y[(N_x + 1) * i:(N_x + 1) * (i + 1)], z[(N_x + 1) * i:(N_x + 1) * (i + 1)],
                  'blue')
    for i in range(N_x + 1):
        ax.plot3D(x2[(N_y + 1) * i:(N_y + 1) * (i + 1)], y2[(N_y + 1) * i:(N_y + 1) * (i + 1)],
                  z2[(N_y + 1) * i:(N_y + 1) * (i + 1)], 'blue')
    plt.show()


def create_u(u):
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
    draw_u(u)


create_u(u_mesh)
