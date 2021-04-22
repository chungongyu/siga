# -*- coding: utf-8 -*-

import sys
from itertools import ifilter,imap,groupby
import string

from fasta import *
from dna_utils import reverse_complement

def contig_read(f, fmt='txt'):
    if fmt == 'txt':
        for contig in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
            yield 'dummy', contig
    else:
        for name, contig in fasta_read(f):
            yield name, contig

def contig_find(ref_tbl, mapping_tbl, contig):
    find = False
    x = 0
    for i, ref in enumerate(ref_tbl):
        r = ref.find(contig, x)
        while r != -1:
            for k in xrange(len(contig)):
                mapping_tbl[i][r+k] = 1
            find = True
            r = ref.find(contig, r + 1)
    return find

if __name__ == '__main__':
    # load reference
    ref = []
    with file(sys.argv[2], 'r') as f:
        print 'ref:'
        for name, seq in fasta_read(f):
            print '%s: %d' % (name, len(seq))
            ref.append(seq)
            #break
    fmt = 'txt'
    if len(sys.argv) > 3:
        fmt = sys.argv[3]
    if len(sys.argv) > 4:
        unmatched_file = open(sys.argv[4],'w')
    else:
        unmatched_file = open('unmatched_contigs','w')

    contig_number = 0
    unnatched_contig = 0
    matched_contig = 0
    contig_total_length = 0
    length_list = []
    mapping_tbl = [[0]*len(i) for i in ref]
    for name, contig in contig_read(sys.stdin, fmt):
        contig_number += 1
        contig_total_length += len(contig)
        length_list.append(len(contig))
        #if contig_find(ref, mapping_tbl, contig) or contig_find(ref, mapping_tbl, reverse_complement(contig)):
        #    matched_contig += 1
        #else:
        #    unmatched_file.write('>' + name + '\n')
        #    unmatched_file.write(contig + '\n')
    print "contig_number: %d"%contig_number
    print "matched_contig: %d"%matched_contig
    print "unmatched_contig: %d"%(contig_number - matched_contig)
    length_list.sort(reverse=True)
    length_sum = 0
    flag = True
    for i in xrange(len(length_list)):
        length_sum += length_list[i]
        if(flag and length_sum > contig_total_length / 2):
            print "N50: %d"%length_list[i]
            flag = False
        if(length_sum > contig_total_length * 9 / 10):
            print "N90: %d"%length_list[i]
            break
    print "MAX_contig: %d"%length_list[0]
    if len(sys.argv) <= 4 or sys.argv[4] != 'disable-mapping':
        for i, mapping in enumerate(mapping_tbl):
            k = 0
            for p, items in groupby(mapping, lambda x: x):
                l = len(list(items))
                print '%d\t%d\t%d\t%d\t%d' % (i, p, k, k + l, l)
                k += l

    unmatched_file.close()
