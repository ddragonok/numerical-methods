import sympy
import numpy
import matplotlib.pyplot as plt

x = sympy.symbols('x')
inp_x = [0, 1, 2, 3, 4]
inp_y = [1, 0.86603, 0.5, 0, -0.5]
n = len(inp_x)

h = numpy.zeros(n - 1)
alpha = numpy.zeros(n - 1)
b = numpy.zeros(n - 1)
d = numpy.zeros(n - 1)
c = numpy.zeros(n)

#  Coefficients of the Tridiagonal Matrix Algorithm
u = numpy.ones(n)
m = numpy.zeros(n)
z = numpy.zeros(n)


def find_coef(inp_x, inp_y):
    for i in range(n - 1):
        h[i] = inp_x[i + 1] - inp_x[i]
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (inp_y[i + 1] - inp_y[i]) - (3 / h[i - 1]) * (inp_y[i] - inp_y[i - 1])
    for i in range(1, n - 1):
        u[i] = 2 * (inp_x[i + 1] - inp_x[i - 1]) - h[i - 1] * m[i - 1]
        m[i] = h[i] / u[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / u[i]
    for i in range(n - 2, -1, -1):
        c[i] = z[i] - m[i] * c[i + 1]
        b[i] = ((inp_y[i + 1] - inp_y[i]) / h[i]) - (h[i] * (c[i + 1] + 2 * c[i]) / 3)
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])

    f = []
    for i in range(n - 1):
        f.append(inp_y[i] + b[i] * (x - inp_x[i]) + c[i] * (pow(x - inp_x[i], 2)) + d[i] * (pow(x - inp_x[i], 3)))
    print(f)
    return f


def cubic_spline(x_in, f, inp_x):
    maxind = n
    for i in range(n - 1, -1, -1):
        if x_in < inp_x[i]:
            maxind = i
    return f[maxind - 1].subs(x, x_in)


def plot_spline(f, inp_x):
    x_point = []
    y_point = []
    count = 100
    delta = (inp_x[n - 1] - inp_x[0]) / count
    for i in range(count):
        x_point.append(inp_x[0] + i*delta)
        y_point.append(cubic_spline(inp_x[0] + i*delta, f, inp_x))
    plt.title('Cubic Spline')
    plt.grid(True)
    plt.plot(x_point, y_point)
    plt.show()


function = find_coef(inp_x, inp_y)
plot_spline(function, inp_x)
