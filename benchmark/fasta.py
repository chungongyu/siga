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
    if len(seq) > 0:
        yield name, make_seq(seq)

if __name__ == '__main__':
    import sys

    def make_test(*args):
        for name, seq in fasta_read(sys.stdin):
            print name
            print seq

    def make_file(*args):
        if len(args) < 2:
            sys.exit(1)

        read1, read2 = file(args[0], 'wb'), file(args[1], 'wb')
        for number, line in enumerate(ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin))):
            if number % 2 == 0:
                print>>read1, '>R%d' % (number/2)
                print>>read1, line
            else:
                print>>read2, '>R%d' % (number/2)
                print>>read2, line

    cmd_tbl = {'make_test': make_test, 'make_file': make_file}
    if len(sys.argv) < 2 or sys.argv[1] not in cmd_tbl:
        print '%s <%s>' % (sys.argv[0], '|'.join(cmd_tbl.keys()))
        sys.exit(1)
    cmd_tbl[sys.argv[1]](*sys.argv[2:])
