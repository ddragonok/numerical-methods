import sympy as sym

a = -1
b = 1
h = 0.5
n = int(abs(a - b) / h)
inp_x = sym.Symbol('x')
inp_y = inp_x / (2*inp_x + 5)


def left_rectangle(x, y, start, step, count):
    f = 0
    for i in range(1, count + 1):
        f += y.subs(x, start + (i-1)*step) * step
    print(f)


def right_rectangle(x, y, start, step, count):
    f = 0
    for i in range(1, count + 1):
        f += y.subs(x, start + i*step) * step
    print(f)


def middle_rectangle(x, y, start, step, count):
    f = 0
    for i in range(1, count + 1):
        f += y.subs(x, ((start + i*step) + (start + (i-1)*step)) / 2) * step
    print(f)


def trapezium(x, y, start, end, step, count):
    f = 0
    summa = 0
    for i in range(1, count):
        summa += y.subs(x, start + i*step)
    f += (summa + (y.subs(x, start) + y.subs(x, end)) / 2) * step
    print(f)


def simpson(x, y, start, end, step, count):
    f = 0
    sum1 = 0
    sum2 = 0
    for i in range(1, count):
        sum1 += y.subs(x, start + i * step)
        sum2 += y.subs(x, ((start + i*step) + (start + (i-1)*step)) / 2)
    f += (2*sum1 + 4*sum2 + y.subs(x, start) + y.subs(x, end)) * (step / 6)
    print(f)


left_rectangle(inp_x, inp_y, a, h, n)
right_rectangle(inp_x, inp_y, a, h, n)
middle_rectangle(inp_x, inp_y, a, h, n)
trapezium(inp_x, inp_y, a, b, h, n)
simpson(inp_x, inp_y, a, b, h, n)
