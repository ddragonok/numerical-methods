import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
y = sym.Symbol('y')
z = sym.Symbol('z')
f = (3*y + z*(x - 2)) / (x-2)**2
a = 3
b = 4
y1 = 2
z1 = 2
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def runge_kutta(h, a, b, y1, z1):
    x_array = [a]
    y_array = [y1]
    z_array = [z1]
    v = int((b - a) / h)
    for i in range(1, v + 1):
        k1 = h * z_array[-1]
        l1 = h * f.subs([(x, x_array[-1]), (y, y_array[-1]), (z, z_array[-1])])
        k2 = h * (z_array[-1] + l1 / 2)
        l2 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k1 / 2), (z, z_array[-1] + l1 / 2)])
        k3 = h * (z_array[-1] + l2 / 2)
        l3 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k2 / 2), (z, z_array[-1] + l2 / 2)])
        k4 = h * (z_array[-1] + l3)
        l4 = h * f.subs([(x, x_array[-1] + h), (y, y_array[-1] + k3), (z, z_array[-1] + l3)])
        delta_y = 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        delta_z = 1 / 6 * (l1 + 2 * l2 + 2 * l3 + l4)
        x_array.append(x_array[-1] + h)
        y_array.append(y_array[-1] + delta_y)
        z_array.append(z_array[-1] + delta_z)
    return x_array, y_array


runge1 = runge_kutta(h1, a, b, y1, z1)
runge2 = runge_kutta(h2, a, b, y1, z1)
draw_graph(runge1, runge2)
