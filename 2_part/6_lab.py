import numpy as np
import math

a = np.array(
    [[5, -4, 7],
     [-4, -3, 4],
     [7, 4, -1]]
    )
N = len(a[0])


def output(A_next, H_comp, n):
    for i in range(n):
        print('l', i, ' = ', A_next[i][i], sep='')
    for j in range(n):
        for i in range(n):
            H_comp[i][j] /= H_comp[n - 1][j]
    for i in range(n):
        print('vector', i, ' = ', np.transpose(H_comp)[i],  sep='')


def triangle_maximum(A, n):
    maximum = 0
    k = 0
    ind1, ind2 = 0, 0
    for i in range(n):
        for j in range(1 + k, n):
            if abs(A[i][j]) > maximum:
                maximum = abs(A[i][j])
                ind1, ind2 = i, j
        k += 1
    return maximum, ind1, ind2


def jacobi(A, n):
    eps = 0.0001
    maximum, ind1, ind2 = triangle_maximum(A, n)

    H = np.zeros((n, n))
    H_comp = np.zeros((n, n))
    for i in range(n):
        H[i][i] = 1
        H_comp[i][i] = 1

    while maximum > eps:
        fi = (math.atan((2 * A[ind1][ind2]) / (A[ind1][ind1] - A[ind2][ind2]))) / 2
        H[ind1][ind1] = math.cos(fi)
        H[ind2][ind2] = math.cos(fi)
        H[ind1][ind2] = -math.sin(fi)
        H[ind2][ind1] = math.sin(fi)

        A_next = np.dot(np.dot(np.transpose(H), A), H)
        H_comp = np.dot(H_comp, H)

        maximum, ind1, ind2 = triangle_maximum(A, n)

        A = np.copy(A_next)

    output(A_next, H_comp, n)


jacobi(a, N)
