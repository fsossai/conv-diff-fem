import numpy as np

filename = '../inputs/mat63.txt'

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
y = read_txt_vec('../vec_y.txt', A.shape[0])
