import matplotlib.pylab as plt
import scipy.sparse as sparse
import numpy as np
import sys

if len(sys.argv) < 2:
    print('Please specify the file name of a sparse CSR matrix.')
    sys.exit(-1)

filename = sys.argv[1]
coef = []
row_ind = []
col_ind = []

print('Reading file...')
with open(filename, 'r') as f:
    s = f.readline().split()
    N = int(s[0])       # Number of rows
    NZ = int(s[1])      # Number of non-zero entries
    
    for line in f:
        s = line.split()
        row_ind.append(int(s[0]) - 1)
        col_ind.append(int(s[1]) - 1)
        coef.append(float(s[2]))

assert NZ == len(coef)

print('Building CSR matrix...')
A = sparse.csr_matrix((coef, (row_ind, col_ind)), shape=(N,N))

print('Creating the plot...')
plt.spy(A)
plt.xticks([])
plt.yticks([])
#plt.show()
outname = filename + '.png'
plt.savefig(outname)

print('Done')
print(f'Plot saved to {outname}')