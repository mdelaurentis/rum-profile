import matplotlib
matplotlib.use('pdf')

import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import os
import argparse

def read_line(f):
    line = f.readline()
    if line is None:
        return (None, None)
    else:
        return (int(x) for x in line.split())

def seqnum_to_alns_iter(filename):
    with open(filename, 'r') as fh:
        for line in fh:
            (alns, seqnum) = line.split()
            yield (seqnum, alns)

def seqnum_to_old_new_iter(old_iter, new_iter):
    
    old_rec = old_iter.next()
    new_rec = new_iter.next()

    while old_rec is not None and new_rec is not None:
        (old_seqnum, old_alns) = old_rec
        (new_seqnum, new_alns) = old_rec
        seqnum = min(old_seqnum, new_seqnum)
        result = (seqnum, 0, 0)

        if (old_seqnum == seqnum):
            result[1] = old_alns
            old_rec = old_iter.next()
        if (new_seqnum == seqnum):
            result[1] = new_alns
            new_rec = new_iter.next()
        yield result

    for rec in old_iter:
        yield (rec[0], rec[1], 0)

    for rec in new_iter:
        yield (rec[0], 0, rec[1])

def mapping_count_iter(filename):
    exhausted = False

    recs = rec_iter(filename)

    with open(filename, 'r') as fh:
        last_rec = recs.next()

        while seqnum is not None:
            counter += 1

            if seqnum == counter:
                result = (alns, seqnum)
                (alns, seqnum) = read_line(f)
            else:
                result = (0, counter)
                yield (0, counter)

            print "Returning " + str(result)
            yield result

colors = ['r', 'b']

def acc_changes(old_filename, new_filename):

    print "Reading %s and %s" % (old_filename, new_filename)

    old_counts = mapping_count_iter(old_filename)
    new_counts = mapping_count_iter(new_filename)

    changes = {}

    while (True):
        old_rec = old_counts.next()
        new_rec = new_counts.next()
        if old_rec is None or new_rec is None:
            break
        (old_alns, old_seqnum) = old_rec
        (new_alns, new_seqnum) = new_rec
        if ((old_seqnum % 100000) == 0):
            print str(old_seqnum)
#        if old_seqnum == 10000000:
#            break
        if (old_seqnum != new_seqnum):
            raise Exception("Different sequence numbers")
        if old_alns not in changes:
            changes[old_alns] = {}
        if new_alns not in changes[old_alns]:
            changes[old_alns][new_alns] = 0
        changes[old_alns][new_alns] += 1

    table = []

    for old_count in sorted(changes):
        for new_count in sorted(changes[old_count]):
            table.append((old_count, new_count, changes[old_count][new_count]))
    return np.array(table)

table_filename='table.npy'

def main():

    if not os.path.isfile(table_filename):
        table = acc_changes(sys.argv[1], sys.argv[2])
        print table
        np.save(table_filename, table)

    table = np.load(table_filename)
    
    x = table[:,0]
    y = table[:,1]
    C = table[:,2]
    plt.hexbin(x, y, C=C, bins='log',gridsize=50)
    plt.colorbar()
    plt.xlabel('alignments in v2.0.2')
    plt.ylabel('alignments in v2.0.3')
    plt.savefig('heatmap')



main()
