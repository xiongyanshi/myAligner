#! /usr/bin/env python

"""
Smith-Waterman pair-wise alignment.
by Yanshi Xiong.
"""

import sys
import argparse
import numpy as np
from lib.seqio import Fasta
from lib.readmat import readmat
from lib.readgap import readgap


def tailn(l, x):
    """
    return the number of continuous element x from the right side of l
    tail([0,1,2,3,4,3,3,3,3], 3) = 4
    """
    count = 0
    for i in l[::-1]:
        if i != x:
            break
        count += 1
    return count

class Align:
    def __init__(self, seq1, seq2, submat, gap_open, gap_extend, algo):
        self.seq1 = '-' + seq1.upper()
        self.seq2 = '-' + seq2.upper()
        self.submat = submat
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.algo = algo

        Align.__buildmat(self)
        Align.__traceback(self)
        Align.__makeprintable(self)

    def __buildmat(self):
        self.scoremat = np.empty([len(self.seq1), len(self.seq2)],
                                 dtype=np.int32)
        self.tracemat = np.empty([len(self.seq1), len(self.seq2)],
                                 dtype=np.int8)

        for i, p in enumerate(self.seq1):
            for j, q in enumerate(self.seq2):
                if i == 0:
                    if self.algo == 'sw':
                        self.scoremat[i,j] = 0
                    elif self.algo == 'nw':
                        self.scoremat[i,j] = j*self.gap_extend + self.gap_open
                    self.tracemat[i,j] = 1
                    continue
                if j == 0:
                    if self.algo == 'sw':
                        self.scoremat[i,j] = 0
                    elif self.algo == 'nw':
                        self.scoremat[i,j] = i*self.gap_extend + self.gap_open
                    self.tracemat[i,j] = 2
                    continue

                # from diagnoal
                d = self.scoremat[i-1, j-1] + self.submat[p][q]

                # from left side
                lgap_len = tailn(self.tracemat[i,:j], 1) + 1   # +1 for itself
                lgap_penalty = lgap_len * self.gap_extend + self.gap_open
                l = self.scoremat[i, j-1] + lgap_penalty

                # from up side
                ugap_len = tailn(self.tracemat[:i,j], 2) + 1   # +1 for itself
                ugap_penalty = ugap_len * self.gap_extend + self.gap_open
                u = self.scoremat[i-1, j] + ugap_penalty

                if self.algo == 'sw':
                    self.scoremat[i, j] = max(d,l,u,0)
                    self.tracemat[i, j] = [d, l, u, 0].index(max(d,l,u,0))
                elif self.algo == 'nw':
                    self.scoremat[i, j] = max(d,l,u)
                    self.tracemat[i, j] = [d, l, u].index(max(d,l,u))

    def __traceback(self):
        self.maxscore = self.scoremat.max()
        self.maxscorei = np.unravel_index(self.scoremat.argmax(),
                                               self.scoremat.shape)
        def reachEnd(self, i, j):
            if self.algo == 'sw':
                return self.scoremat[i,j] == 0     # 0-scored cell is end.
            else:
                return (i,j) == (0,0)              # [0, 0] cell is end.

        if self.algo == 'sw':
            self.end = self.maxscorei              # start loc of tracing back
        elif self.algo == 'nw':
            self.end = len(self.seq1)-1, len(self.seq2)-1

        self.pathcode = ''
        i, j = self.end
        while not reachEnd(self, i, j):
            #print('[%s,%s]<' % (i, j), end='')
            direction = str(self.tracemat[i,j])
            self.pathcode = direction + self.pathcode
            if direction == '0':
                i -= 1
                j -= 1
            elif direction == '1':
                j -= 1
            elif direction == '2':
                i -= 1
        self.start = [i,j]
        # the end of tracing back journey,
        # start of alignment segment.

    def __makeprintable(self):
        i, j = self.start  # start of aligned base
        i = i+1
        j = j+1

        # unaligned sequence at left side
        leftsize = max(i, j) - 1                   # -1 : cut '-' prefix
        top    = '{:>{leftsize}}'.format(self.seq1[1:i], leftsize=leftsize)
        middle = ' ' * leftsize
        bottom = '{:>{leftsize}}'.format(self.seq2[1:j], leftsize=leftsize)

        # aligned body
        for path in self.pathcode:
            #print('[%s,%s]:%s' % (i, j, path))
            if path == '0':
                top    += self.seq1[i]
                middle += '|' if self.seq1[i] == self.seq2[j] else '*'
                bottom += self.seq2[j]
                i += 1
                j += 1
            if path == '1':
                top    += '-'
                middle += ' '
                bottom += self.seq2[j]
                j += 1
            if path == '2':
                top    += self.seq1[i]
                middle += ' '
                bottom += '-'
                i += 1

        # unaligned sequence at right side
        m, n = self.end
        m += 1
        n += 1
        rightsize = max(len(self.seq1) - m, len(self.seq2) - n)
        top    += '{:<{rightsize}}'.format(self.seq1[m:], rightsize=rightsize)
        middle += ' ' * rightsize
        bottom += '{:<{rightsize}}'.format(self.seq2[n:], rightsize=rightsize)

        top = top.replace(' ', '-')       # for un-aligned part
        bottom = bottom.replace(' ', '-')
        res = "\n  %s\n  %s\n  %s\n" % (top, middle, bottom)
        self.alignedseq1 = top
        self.alignedseq2 = bottom
        self.printable = res

    def print_align(self):
        print(self.printable)

    def print_scoremat(self):
        print(self.scoremat)

    def print_tracemat(self):
        print(self.tracemat)

    def __repr__(self):
        return self.printable


def main():
    # pair-wise alignment with smith-waterman algorithm.

    # substitution matrix and gap penalty are must.
    submat     = readmat('./lib/dna_sub.default.mat')
    gap_open   = readgap('./lib/dna_gap.default.txt')['gap_open']
    gap_extend = readgap('./lib/dna_gap.default.txt')['gap_extend']
    algo = 'sw'

    seq1 = sys.argv[1]
    seq2 = sys.argv[2]

    # ready to go
    align = Align(seq1, seq2, submat, gap_open, gap_extend, algo)
    align.print_align()
    #align.print_scoremat()
    #align.print_tracemat()

    return 0

if __name__ == '__main__':
    main()
