# -*- coding: utf-8 -*-

from itertools import ifilter,imap
import string

def fasta_read(f):
    trans_tbl = {
            'U': 'T', 
            'R': 'A', 
            'Y': 'C', 
            'S': 'G', 
            'W': 'A', 
            'K': 'G', 
            'M': 'A', 
            'B': 'C', 
            'D': 'A', 
            'H': 'C', 
            'V': 'G', 
            'N': 'T', 
            'A': 'A', 
            'C': 'C', 
            'G': 'G', 
            'T': 'T' 
            }

    def make_seq(seq):
        return ''.join(imap(lambda x: trans_tbl[x] if trans_tbl.has_key(x) else x, seq))

    name, seq = '', ''
    state = 0
    for line in ifilter(lambda x: len(x)>0, imap(string.strip, f)):
        if state == 0:
            if line.startswith('>'):
                name = line[1:]
                state = 1
            else:
                raise Exception('invalid data')
        elif state == 1:
            if line.startswith('>'):
                yield name, make_seq(seq)
                name = line[1:]
                seq = ''
            else:
                seq += line
            state = 2
        elif state == 2:
            state = 3
        else:
            state = 0
    if len(seq) > 0:
        yield name, seq

if __name__ == '__main__':
    import sys
    table = {'A', 'T', 'C', 'G'}
    for name, seq in fasta_read(sys.stdin):
        a = seq.count('A')
        t = seq.count('T')
        c = seq.count('C')
        g = seq.count('G')
        if a + t + c + g != len(seq):
            print seq
