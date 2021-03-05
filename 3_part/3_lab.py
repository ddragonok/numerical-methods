import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
y = sym.Symbol('y')
f = ((y**2)*sym.log(x) - y) / x
a = 1
b = 2
y1 = 0.5
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def runge_kutta(h, a, y1):
    x_array = [a]
    y_array = [y1]
    for i in range(1, 4):
        k1 = h * f.subs([(x, x_array[-1]), (y, y_array[-1])])
        k2 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k1 / 2)])
        k3 = h * f.subs([(x, x_array[-1] + h / 2), (y, y_array[-1] + k2 / 2)])
        k4 = h * f.subs([(x, x_array[-1] + h), (y, y_array[-1] + k3)])
        delta = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x_array.append(x_array[-1] + h)
        y_array.append(y_array[-1] + delta)
    return x_array, y_array


def adams(h, x_array, y_array, b):
    while x_array[-1] <= b:
        y_array.append(y_array[-1] + (h / 24)*(55 * f.subs([(x, x_array[-1]), (y, y_array[-1])])
                                               - 59 * f.subs([(x, x_array[-2]), (y, y_array[-2])])
                                               + 37 * f.subs([(x, x_array[-3]), (y, y_array[-3])])
                                               - 9 * f.subs([(x, x_array[-4]), (y, y_array[-4])])))
        x_array.append(x_array[-1] + h)
    return x_array, y_array


ad1 = runge_kutta(h1, a, y1)
ad1 = adams(h1, ad1[0], ad1[1], b)
ad2 = runge_kutta(h2, a, y1)
ad2 = adams(h2, ad2[0], ad2[1], b)
draw_graph(ad1, ad2)
