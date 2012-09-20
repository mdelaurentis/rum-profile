import matplotlib
matplotlib.use('pdf')

import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import os
import argparse
import sqlite3
import itertools

table = []

with open('out') as f:
    f.next()
    for line in f:
        row = tuple([int(x) for x in line.split()])
        (count, old_alns, new_alns) = row
        table.append(row)

print str(table)
table = np.array(table)
C = table[:,0]
x = table[:,1]
y = table[:,2]
plt.hexbin(x, y, C=C, bins='log', gridsize=50)
#plt.colorbar()
#plt.xlabel('alignments in v2.0.2')
#plt.ylabel('alignments in v2.0.3')
#plt.xlim(0, 1000)
#plt.ylim(0, 1000)
plt.savefig('heatmap')

