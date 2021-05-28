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
single = p * n
ideal = range(1, max(data[0])+1)
speedup = [single / x for x in data[1]]

r = range(1, len(ideal)+1)
plt.plot(data[0], speedup, '.-', label='Benchmarked')
plt.plot(r, ideal, '--', label='Ideal')
plt.xticks(r, r)
plt.title('BiCGSTAB - OpenMP strong scaling')
plt.xlabel('Number of threads')
plt.ylabel('Speedup')
plt.grid(True)
plt.legend()
plt.show()
