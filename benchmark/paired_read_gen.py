# -*- coding: utf-8 -*-

import string
import sys
import random
import numpy

from fasta import *

random.seed()

K = int(sys.argv[2])
coverage = K
if len(sys.argv) > 3:
    coverage = int(sys.argv[3])
insertSize = 1000
if len(sys.argv) > 4:
    insertSize = int(sys.argv[4])
sigma = 0
if len(sys.argv) > 5:
    sigma = int(sys.argv[5])
deltaSet = []
if sigma != 0:
    deltaSet = numpy.random.normal(0, sigma, 10000)

ref = ''
with file(sys.argv[1], 'rb') as f:
    for _, seq in fasta_read(f):

        ref = seq.upper()

        ref_length = len(ref)
        if ref_length < 2 * K + insertSize + 5 * sigma:
            continue
        postPos = set()
        readDict = {}
        for i in xrange(ref_length * coverage / (2 * K)):
            pos = random.randint(0,ref_length - 2 * K - insertSize)
            while pos in postPos:
                pos = random.randint(0,ref_length - 2 * K - insertSize)
            postPos.add(pos)
            delta = 0
            if sigma != 0:
                t = random.randint(0, 9999)
                delta = int(deltaSet[t])
            if pos + 2*K + insertSize + delta > len(ref):
                continue
            r1 = ref[pos:pos + K]
            r2 = ref[pos + K + insertSize + delta: pos + 2 * K + insertSize + delta]
            read = r1 + r2
            if readDict.has_key(read) or read.count('N') >= 10:
                continue
            else:
                readDict[read] = 0
                print r1
                print r2
        '''
        sp = sorted(postPos)
        for i in range(len(sp) - 1):
            if sp[i + 1] - sp[i] > 1000:
                print sp[i],sp[i + 1]
        '''
