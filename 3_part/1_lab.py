import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
y = sym.Symbol('y')
f = x*x*(y*y + 1)
a = 0
b = 1
y0 = 0
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def euler(h, a, b, y0):
    x_array = [a]
    y_array = [y0]
    v = int((b - a) / h)
    for i in range(1, v + 1):
        x_array.append(x_array[-1] + h)
        y_array.append(y_array[-1] + h*f.subs([(x, x_array[-1]), (y, y_array[-1])]))
    return x_array, y_array


def upgrade_euler(h, a, b, y0):
    x_array = [a]
    y_array = [y0]
    v = int((b - a) / h)
    for i in range(1, v + 1):
        arg = y_array[-1] + h*f.subs([(x, x_array[-1]), (y, y_array[-1])])
        x_array.append(x_array[-1] + h)
        mult = f.subs([(x, x_array[-2]), (y, y_array[-1])]) + f.subs([(x, x_array[-1]), (y, arg)])
        y_array.append(y_array[-1] + (h / 2)*mult)
    return x_array, y_array


eul1 = euler(h1, a, b, y0)
eul2 = euler(h2, a, b, y0)
draw_graph(eul1, eul2)

up_eul1 = upgrade_euler(h1, a, b, y0)
up_eul2 = upgrade_euler(h2, a, b, y0)
draw_graph(up_eul1, up_eul2)



