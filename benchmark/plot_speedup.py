import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import scipy.optimize

def amdahls_law(x, s):
    return 1 / (s + (1 - s)/x)

def modified_amdahls_law(x, s, v):
    return 1 / (s + (1 - s - v)/x + v * x)

if len(sys.argv) < 2:
    print('ERROR: missing file name')
    sys.exit(-1)

filename = sys.argv[1]
data = pd.read_csv(filename, sep='\t', header=None)

# Considering only the first 'sys.argv[2]' rows if specified
if (len(sys.argv) == 3):
    data = data[0:int(sys.argv[2])]
    
_, (n,p) = next(data.iterrows())
single = p * n
ideal = range(1, max(data[0])+1)
speedup = [single / x for x in data[1]]

# calculating the fraction of serial program according
# to Amdhal's Law
par, _ = scipy.optimize.curve_fit(amdahls_law, data[0], speedup)
#par, _ = scipy.optimize.curve_fit(modified_amdahls_law, data[0], speedup)
s = par[0]
#v = par[1]
print('s =', s)
#print('v =', v)

r = range(1, len(ideal)+1)
plt.plot(data[0], speedup, '.-', label='Benchmarked')
#plt.plot(r, ideal, '--', label='Ideal')
plt.plot(r, amdahls_law(np.array(r), par[0]), '--', label=f"Amdahl's Law, s={s:.3}")
#plt.plot(r, modified_amdahls_law(np.array(r), par[0]), '--', label=f"Modified Amdahl's Law, s={s:.3}")
plt.xticks(r, r)
plt.title('BiCGSTAB - OpenMP strong scaling')
plt.xlabel('Number of threads')
plt.ylabel('Speedup')
plt.grid(True)
plt.legend()
plt.show()
