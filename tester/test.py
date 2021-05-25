"""
Simple utilities to parse sparse matrix and vectors.
The purpose is to asses whether some performed calculations
are correct.
"""

import numpy as np
from scipy.sparse.linalg import bicgstab as scipy_bicgstab
import matplotlib.pyplot as plt
from time import time

def bicgstab(A, b, toll=1e-5, max_it=100):
    N, _ = A.shape
    x = [np.zeros((N,))]
    r = [b - A.dot(x[0])]
    p = r[0].copy()
    #rs = np.ones((N,))
    rs = r[0].copy()

    j = 0
    while np.linalg.norm(r[j]) > toll and j < max_it:
        Ap = A.dot(p)
        alpha = r[j].dot(rs) / Ap.dot(rs)
        s = r[j] - alpha*Ap
        As = A.dot(s)
        omega = (As.dot(s)) / (As.dot(As))
        x.append(x[j] + alpha*p + omega*s)
        r.append(s - omega*As)
        beta = r[j+1].dot(rs) / r[j].dot(rs) * alpha / omega
        p = r[j+1] + beta*(p - omega*Ap)
        j += 1
    return x,r,j

def solve(A, b, toll=1e-05, max_it=100, plot=True):
    t = time()
    x, r, it = bicgstab(A, b, toll, max_it)
    t = time() - t
    print('Elapsed time (raw bicgstab):', t)

    # solving with matrix inversion
    t = time()
    sol = np.linalg.inv(A).dot(b)
    t = time() - t
    print('Elapsed time (matrix inv):', t)

    print('Error:', np.linalg.norm(b - A.dot(x[-1])))
    print('Residual:', np.linalg.norm(r[-1]))
    print('Min residual:', min([np.linalg.norm(v) for v in r]))
    print('Iteration:', it)
    print()
    print('Solution:', sol)
    print('Found:', x[-1])
    if plot:
        plt.plot(range(len(r)), [np.linalg.norm(v) for v in r], '.-')
        plt.yscale('log')
        #plt.xticks(range(0,it+1,5), range(0,it+1,5))
        plt.grid(True)
        plt.show()

def read_txt_mat(filename):
    """reads a sparse matrix from a text file"""

    with open(filename, 'r') as f:
        header = f.readline()
        N = int(header.split()[0])
        nonzeros = int(header.split()[1])
        A = np.zeros((N,N))
        for line in f:
            s = line.split()
            i = (int(s[0]) - 1, int(s[1]) - 1)
            A[i] = float(s[2])
    return A

def read_txt_vec(filename, N):
    """reads a sparse vector from a text file"""

    with open(filename, 'r') as f:
        v = np.zeros((N,))
        for line in f:
            s = line.split()
            v[int(s[0]) - 1] = float(s[1])
    return v

A = read_txt_mat('../inputs/mat63.txt')
N = A.shape[0]
#y = read_txt_vec('../solution.txt', N)
b = np.ones((N,))

t = time()
M = np.eye(N)
#M = np.diag(1 / np.diag(A))
#M = np.linalg.inv(A)
x, _, = scipy_bicgstab(A, b, tol=1e-6, M=M)
t = time() - t
print('scipy:')
print(x)
print(f'Elapsed time:', t)

t = time()
solve(M.dot(A), M.dot(b), toll=1e-6, max_it=N*10, plot=True)
t = time() - t
print(f'Elapsed time:', t)

