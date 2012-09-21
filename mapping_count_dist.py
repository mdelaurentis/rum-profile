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

with open('data/counts') as f:
    f.next()
    for line in f:
        row = tuple([float(x) for x in line.split()])
        (count, old_alns, new_alns) = row
        table.append(row)


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


classes = ['none', 'unique', 'non-unique']

with open('data/classes') as f:
    f.next()
    table = []

    for line in f:
        row = tuple([float(x) for x in line.split()])
        table.append(row)

    table = np.array(table, dtype=[
            ('reads', int),
            ('old_class_num', int),
            ('new_class_num', int)])

    total = np.sum(table['reads'])
    
    table = np.sort(table, order=('reads'))
    table = table[::-1]

    acc = 0.0

    num_changed = 0.0

    print "Total is " + str(total)

    for row in table:
        reads     = row['reads']
        old_class = classes[row['old_class_num']]
        new_class = classes[row['new_class_num']]
        if old_class != new_class:
            num_changed += reads
        pct = reads * 100.0 / total
        acc += pct
        print "%s to %s: %d (%.5f%%) (%.5f%%)" % (old_class, new_class, reads, pct, acc)
    
    print "%f: %f" % (num_changed, 100.0 * (num_changed / total))
    num_same = total - num_changed
    print "%f: %f" % (num_same, 100.0 * (num_same / total))

    print table
    
