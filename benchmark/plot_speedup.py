import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import scipy.optimize

def amdahls_law(x, s):
    return 1 / (s + (1 - s)/x)

if len(sys.argv) < 2:
    print('ERROR: missing file name')
    sys.exit(-1)

filename = sys.argv[1]
data = pd.read_csv(filename, sep='\t', header=None)
_, (n,p) = next(data.iterrows())
single = p * n
ideal = range(1, max(data[0])+1)
speedup = [single / x for x in data[1]]

# calculating the fraction of serial program according
# to Amdhal's Law
s, _ = scipy.optimize.curve_fit(amdahls_law, data[0], speedup)
s = s[0]
print('s =', s)

r = range(1, len(ideal)+1)
plt.plot(data[0], speedup, '.-', label='Benchmarked')
plt.plot(r, ideal, '--', label='Ideal')
plt.plot(r, amdahls_law(np.array(r), s), '--', label=f"Amdahl's Law, s={s:.3}")
plt.xticks(r, r)
plt.title('BiCGSTAB - OpenMP strong scaling')
plt.xlabel('Number of threads')
plt.ylabel('Speedup')
plt.grid(True)
plt.legend()
plt.show()
