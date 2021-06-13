import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 3:
    print('Please specify an input file name AND the number of nodes.')
    sys.exit(-1)

filename = sys.argv[1]
N = int(sys.argv[2])
conn = [set() for _ in range(N)]

with open(filename, 'r') as f:
    f.readline()
    for line in f:
        s = line.split()
        a, b, c = int(s[1]) - 1, int(s[2]) - 1, int(s[3]) - 1

        conn[a].add(b); conn[a].add(c)
        conn[b].add(a); conn[b].add(c)
        conn[c].add(a); conn[c].add(b)

deg = list(map(len, conn))
avg = np.array(deg).mean()
print('Average node degree:', avg)

plt.bar(*np.unique(deg, return_counts=True))
plt.xticks()
plt.title(f'Degree distribution of the topology matrix\n{filename}')
plt.xlabel('Number of connections (degree)')
plt.ylabel('Count')
plt.show()
