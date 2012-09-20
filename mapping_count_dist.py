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

def seqnum_iter(filename):
    pat = re.compile("^seq\.(\d+)(a|b)?")

    last_seqnum    = None
    last_direction = None

    with open(filename) as f:
        for line in f:
            m = pat.match(line)
            if m is None:
                raise Exception("Bad line: " + line)
            seqnum    = int(m.group(1))
            direction = m.group(2)
            if (seqnum         == last_seqnum and
                direction      == 'b'         and
                last_direction == 'a'):
                pass
            else:
                yield seqnum
            last_seqnum    = seqnum
            last_direction = direction

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

colors = ['r', 'b']

def acc_changes(old_filename, new_filename):

    recs = seqnum_to_old_new_iter(
        seqnum_to_alns_iter(old_filename),
        seqnum_to_alns_iter(new_filename))

    changes = {}
    counter = 0
    for rec in recs:
        counter += 1
        (seqnum, old_alns, new_alns) = rec
        if ((counter % 100000) == 0):
            print str(old_seqnum)
        if old_seqnum == 100000:
            break
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

db_filename = 'aln_count_delta.db'

def load_seqnums(conn, directory, label):
    try:
        conn.execute('''CREATE TABLE %s_seqnum_alns (seqnum integer)''' % label)
    except:
        print "Looks like we've already loaded seqnums from %s" % directory
        return

    old_unique = seqnum_iter(directory + "/RUM_Unique")
    old_nu     = seqnum_iter(directory + "/RUM_NU")

    print "  Reading RUM_Unique and RUM_NU"
    c = conn.cursor()
    for seqnum in itertools.chain(old_unique, old_nu):
        c.execute('INSERT INTO %s_seqnum_alns VALUES(?)' % label, (seqnum,))

    c.execute('SELECT count(1) FROM %s_seqnum_alns' % label)
    row = c.fetchone()
    print "    Loaded %d alignments" % row[0]
    c.close()

def count_alns(conn, label):
    try:
        conn.execute('''CREATE TABLE %s_seqnum_to_aln_count (seqnum integer primary key, alns integer)''' % label)
    except:
        print '  Looks like we\'ve already counted alignments for %s' % label
        return

    c1 = conn.cursor()
    c2 = conn.cursor()
    print "  Counting alignments for each read"
    for row in c1.execute('''SELECT   seqnum, 
                                     count(seqnum) 
                              FROM   %s_seqnum_alns
                            GROUP BY seqnum''' % label):
        c2.execute('INSERT INTO %s_seqnum_to_aln_count VALUES(?, ?)' % label, row)

    c2.execute('SELECT count(1) FROM %s_seqnum_to_aln_count' % label)
    row = c2.fetchone()
    print "    Counted %d distinct reads" % row[0]
    conn.commit()
    c1.close()
    c2.close()

def find_diffs(conn, old_label, new_label):
    c1 = conn.cursor()
    print 'Getting diffs from %s and %s' % (old_label, new_label)

    try:
        c1.execute('CREATE TABLE diffs (reads integer, old_alns integer, new_alns integer)')
    except:
        print '  Looks like we\'ve already found diffs'
        return


    queries = [
        '''    SELECT count(old.seqnum), 
                      old.alns, 
                      new.alns
                 FROM %s_seqnum_to_aln_count old INNER JOIN
                      %s_seqnum_to_aln_count new ON old.seqnum = new.seqnum
             GROUP BY old.alns, new.alns''' % (old_label, new_label),
    
        '''    SELECT count(old.seqnum), 
                      old.alns, 
                      0
                 FROM %s_seqnum_to_aln_count old LEFT JOIN
                      %s_seqnum_to_aln_count new ON old.seqnum = new.seqnum
                WHERE new.seqnum IS NULL
             GROUP BY old.alns, new.alns''' % (old_label, new_label),

        '''    SELECT count(new.seqnum), 
                      0,
                      new.alns
                 FROM %s_seqnum_to_aln_count new LEFT JOIN
                      %s_seqnum_to_aln_count old ON old.seqnum = new.seqnum
                WHERE old.seqnum IS NULL
             GROUP BY old.alns, new.alns''' % (new_label, old_label)
        ]
    
    c2 = conn.cursor()

    for q in queries:
        for row in c1.execute(q):
            c2.execute('INSERT INTO diffs VALUES(?, ?, ?)', row)
    conn.commit()

def load_counts(conn, directory, label):
    print "Loading counts for %s from %s" % (label, directory)
    load_seqnums(conn, directory, label)
    count_alns(conn, label)


def main():

    conn = sqlite3.connect(db_filename)

    old_dir = sys.argv[1]
    new_dir = sys.argv[2]

    load_counts(conn, old_dir, 'old')
    load_counts(conn, new_dir, 'new')
    find_diffs(conn, 'old', 'new')



#    if not os.path.isfile(table_filename):
#        table = acc_changes(sys.argv[1], sys.argv[2])
#        print table
#        np.save(table_filename, table)

#    table = np.load(table_filename)
    
#    x = table[:,0]
#    y = table[:,1]
#    C = table[:,2]
#    plt.hexbin(x, y, C=C, bins='log',gridsize=50)
#    plt.colorbar()
#    plt.xlabel('alignments in v2.0.2')
#    plt.ylabel('alignments in v2.0.3')
#    plt.savefig('heatmap')



main()
