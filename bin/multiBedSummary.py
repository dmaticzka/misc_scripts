#! /usr/bin/env python

from __future__ import absolute_import
from pybedtools import BedTool
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description="""

Given two or mote bed files, ``multiBedSumary.py`` computes the sum of overlapping intervals in every genomic region. The default output of ``multiBedSumary.py`` (a compressed numpy array, .npz) can be used from various tools of the deepseq2 package such as ``plotCorrelation`` or ``plotPCA`` for visualization and diagnostic purposes.

""")

# required = parser.add_argument_group('Required arguments')

parser.add_argument('--regions', '-r',
                    help='BED file containing all regions that should be considered.',
                    required=True)

parser.add_argument('--bedfiles', '-b',
                    help='List of bed files, separated by spaces.',
                    nargs='+',
                    required=True)

parser.add_argument('--outFileName', '-out',
                    help='File name to save the compressed matrix file (npz format)'
                    'needed by the "plotHeatmap" and "plotProfile" tools.',
                    required=True)

parser.add_argument('--labels', '-l',
                    help='User defined labels instead of default labels from '
                    'file names. '
                    'Multiple labels have to be separated by spaces, e.g., '
                    '--labels sample1 sample2 sample3',
                    nargs='+')

args = parser.parse_args()


# count all overlaps with reference intervals
annot = BedTool().annotate(i=args.regions,
                           files=args.bedfiles,
                           counts=True)

# load counts to numpy array
counts = np.loadtxt(annot.fn, usecols=list(range(3, 3 + len(args.bedfiles))))

# combine matrix and lavels, save as npz file
if args.labels is None:
    labels = args.bedfiles
else:
    labels = args.labels
np.savez(args.outFileName, labels=np.array(labels), matrix=counts)
