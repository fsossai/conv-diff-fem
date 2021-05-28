import numpy as np
import scipy.sparse
from scipy.sparse.linalg import bicgstab as scipy_bicgstab

def export(A, filename):
    N = A.shape[0]
    with open(filename, 'w') as f:
        f.write(f'{N} {(A != 0).sum()}\n')
        for i in range(N):
            for j in range(N):
                if A[i,j] != 0:
                    f.write(f'{i+1} {j+1} {A[i,j]}\n')
                

N = 1500
A = scipy.sparse.random(N,N,density=0.01,dtype=np.float64).toarray()
A = A*A
A = A.dot(A.T)
b = np.ones((N,), dtype=np.float64)
print(f'Sparseness: {(A==0).sum()/(N*N)*100:.5} %')
outp = scipy_bicgstab(A,b,tol=1e-6)
print('Computing inversion')
x = outp[0]
x_real = np.linalg.inv(A).dot(b)

print(x[:20])
print()
print(x_real[:20])
print('bicgstab exit code:', outp[1])
export(A, f'mat{N}.txt')
