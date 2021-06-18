import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 2:
    print('Please specify an input file name.')
    sys.exit(-1)

filename = sys.argv[1]

x = []
y = []

with open(filename, 'r') as f:
    N = int(f.readline())
    for line in f:
        s = line.split()
        x.append(float(s[1]))
        y.append(float(s[2]))

x = np.array(x)
y = np.array(y)

# Finding nodes on the boundaries

tol = 1e-5
# Points A and B define a squared
point_A = [x.min(), y.min()]
point_B = [x.max(), y.max()]

I = np.arange(N, dtype=np.int32)
vertical_bound_A = I[np.abs(x-point_A[0]) < tol].tolist()
vertical_bound_B = I[np.abs(x-point_B[0]) < tol].tolist()
horizontal_bound_A = I[np.abs(y-point_A[1]) < tol].tolist()
horizontal_bound_B = I[np.abs(y-point_B[1]) < tol].tolist()
bound = list(set(
    vertical_bound_A + vertical_bound_B +
    horizontal_bound_A + horizontal_bound_B
))
bound = np.array(bound)
bound.sort()

print('Number of boundary nodes:', len(bound))
print('Boundary nodes:')
print(bound)

N = len(x)
fractions = 16
pen_size=0.05
for i in range(fractions):
    i_min = (N // fractions) * i
    i_max = (N // fractions) * (i + 1)
    plt.scatter(x[i_min:i_max], y[i_min:i_max], s=pen_size)
plt.scatter(x[bound], y[bound], s=pen_size)
plt.show()