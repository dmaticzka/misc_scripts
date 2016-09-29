#! /usr/bin/env python

import argparse
# from pybedtools import BedTool

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
                    metavar='sample1 sample2',
                    help='User defined labels instead of default labels from '
                    'file names. '
                    'Multiple labels have to be separated by spaces, e.g., '
                    '--labels sample1 sample2 sample3',
                    nargs='+')

args = parser.parse_args()

# BedTool().annotate
