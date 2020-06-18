# -*- coding: utf-8 -*-

import sys
from itertools import ifilter,imap
import string
import copy

def reverse_complement(read):
    rcread = ''
    for c in read:
        if c == 'A':
            rcread += 'T'
        elif c == 'C':
            rcread += 'G'
        elif c == 'G':
            rcread += 'C'
        elif c == 'T':
            rcread += 'A'
        else:
            rcread += 'A'
    return rcread[::-1]

if __name__ == '__main__':
    fmt = 'txt'
    if len(sys.argv) > 1:
        fmt = sys.argv[1]

    if fmt == 'txt':
        for contig in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
            print reverse_complement(contig)
    elif fmt == 'fasta':
        from fasta import fasta_read
        for name, contig in fasta_read(sys.stdin):
            print '>%s' % (name)
            print reverse_complement(contig)
