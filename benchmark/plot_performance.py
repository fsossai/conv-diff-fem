import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd

if len(sys.argv) < 2:
    print('ERROR: missing file name')
    sys.exit(-1)

filename = sys.argv[1]
data = pd.read_csv(filename, sep='\t', header=None)
_, (n,p) = next(data.iterrows())
single = p / n
ideal = [single * i for i in range(1, max(data[0])+1)]

r = range(1, len(ideal)+1)
plt.plot(data[0], data[1], '.-', label='Benchmarked')
plt.plot(r, ideal, '--', label='Ideal')
plt.xticks(r, r)
plt.title('amxpby - OpenMP strong scaling')
plt.xlabel('Number of threads')
plt.ylabel('Performance [GFlops]')
plt.grid(True)
plt.legend()
plt.show()
