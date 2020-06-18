# -*- coding:utf-8 -*-
import sys
from itertools import ifilter,imap
import string

import graphviz

def lines(f):
    for line in ifilter(lambda x: len(x)>0, imap(string.strip, f)):
        yield line.split('\t')

def main(args):
    color_tbl = {
            '0':'yellow', 
            '1':'gray', 
            '2':'black', 
            '3':'red', 
            '4':'green', 
            '5':'blue'
            }
    red_nodes = set(args.red_node) if args.red_node else set()

    g = graphviz.Digraph(comment='siga')
    for fields in lines(sys.stdin):
        if fields[0] == 'VT':
            attrs = {'style':'filled','color':'red'} if fields[1] in red_nodes else {}
            g.node(fields[1], label='%s:%d' % (fields[1], len(fields[2])), **attrs)
        elif fields[0] == 'ED':
            v1, v2, overlap = fields[1].split(' ', 2)
            xb, xe, xl, yb, ye, yl, rc, _ = map(int, overlap.split())
            attrs = {}
            for i in xrange(2, len(fields)):
                key, typ, val = fields[i].split(':')
                if key == 'CL':
                    attrs['color'] = color_tbl.get(val, 'black')
            if xb == 0 and yb == 0:
                g.edge(v1, v2, label='%d:r:%d' % (xe-xb+1, rc), **attrs)
            elif xb == 0:
                g.edge(v2, v1, label='%d:f:%d' % (xe-xb+1, rc), **attrs)
            elif yb == 0:
                g.edge(v1, v2, label='%d:f:%d' % (xe-xb+1, rc), **attrs)
            else:
                g.edge(v2, v1, label='%d:f:%d' % (xe-xb+1, rc), **attrs)

    print g.source

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser();
    parser.add_argument('-min-overlap', metavar='INT', help='min overlap', type=int)
    parser.add_argument('-red-node', metavar='NODE', help='red nodes', type=str, nargs='+')

    main(parser.parse_args())

