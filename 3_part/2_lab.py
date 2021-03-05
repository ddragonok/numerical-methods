import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
y = sym.Symbol('y')
f = (y / x) + sym.sin(x)*x**2
a = 1
b = 2
y1 = 1.3011687
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def runge_kutta(h, a, b, y1):
    x_array = [a]
    y_array = [y1]
    v = int((b - a) / h)
    for i in range(1, v + 1):
        k1 = h * f.subs([(x, x_array[-1]), (y, y_array[-1])])
        k2 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k1 / 2)])
        k3 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k2 / 2)])
        k4 = h * f.subs([(x, x_array[-1] + h), (y, y_array[-1] + k3)])
        delta = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x_array.append(x_array[-1] + h)
        y_array.append(y_array[-1] + delta)
    return x_array, y_array


runge1 = runge_kutta(h1, a, b, y1)
runge2 = runge_kutta(h2, a, b, y1)
draw_graph(runge1, runge2)