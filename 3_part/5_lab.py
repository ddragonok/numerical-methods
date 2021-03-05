import sympy
import matplotlib.pyplot as plt

x = sympy.symbols('x')
y = sympy.symbols('y')
z = sympy.symbols('z')
f = 2 * y / (x**2 + 1)
a = 0
b = 2
y_a = 0
y_b = 1
h1 = 0.1
h2 = 0.01


def draw_graph(first, second):
    plt.grid(True)
    plt.plot(first[0], first[1], color='red')
    plt.plot(second[0], second[1], color='blue')
    plt.show()


def runge_kutta(h, x_array, y_array, z_array, a, b):
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
    return x_array, y_array, z_array


def find_foo(h, nu, cond_x, cond_z, a, b):
    cond_y = [nu]
    cond_x, cond_y, cond_z = runge_kutta(h, cond_x, cond_y, cond_z, a, b)
    return cond_y[-1], cond_z[-1]


def find_nu(h, nu_1, nu_2, cond_x, cond_z, a, b, y_b):
    eps = 0.001
    array_nu = [nu_1, nu_2]
    y_last, z_last = find_foo(h, nu_1, cond_x, cond_z, a, b)
    foo = [y_last - z_last - y_b]
    if abs(foo[-1]) < eps:
        return nu_1
    y_last, z_last = find_foo(h, nu_2, cond_x, cond_z, a, b)
    foo.append(y_last - z_last - y_b)
    while abs(foo[-1]) > eps:
        array_nu.append(array_nu[-1] - ((array_nu[-1] - array_nu[-2]) * foo[-1]) / (foo[-1] - foo[-2]))
        y_last, z_last = find_foo(h, array_nu[-1], cond_x, cond_z, a, b)
        foo.append(y_last - z_last - y_b)
    return array_nu[-1]


arrayY1 = [find_nu(h1, 1, 2, [a], [y_a], a, b, y_b)]
arrayX1 = [a]
arrayZ1 = [y_a]
runge1 = runge_kutta(h1, arrayX1, arrayY1, arrayZ1, a, b)

arrayY2 = [find_nu(h2, 1, 2, [a], [y_a], a, b, y_b)]
arrayX2 = [a]
arrayZ2 = [y_a]
runge2 = runge_kutta(h2, arrayX2, arrayY2, arrayZ2, a, b)

draw_graph(runge1, runge2)
