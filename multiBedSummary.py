#! /usr/bin/env python

import argparse
#from pybedtools import BedTool

parser = argparse.ArgumentParser(
    description="""

Given two or mote bed files, ``multiBedSumary.py`` computes the sum of overlapping intervals in every genomic region. The default output of ``multiBedSumary.py`` (a compressed numpy array, .npz) can be used from various tools of the deepseq2 package such as ``plotCorrelation`` or ``plotPCA`` for visualization and diagnostic purposes.

""")

required = parser.add_argument_group('Required arguments')

required.add_argument('--regions', '-r',
                      help='BED file containing all regions that should be considered.',
                      nargs='1',
                      required=True)

required.add_argument('--bedfiles', '-b',
                      metavar='FILE1 FILE2',
                      help='List of bed files, separated by spaces.',
                      nargs='+',
                      required=True)

required.add_argument('--outFileName', '-out',
                      help='File name to save the compressed matrix file (npz format)'
                      'needed by the "plotHeatmap" and "plotProfile" tools.',
                      required=True)

optional = parser.add_argument_group('Optional arguments')

optional.add_argument("--help", "-h", action="help",
                      help="show this help message and exit")

# optional.add_argument('--labels', '-l',
#                       metavar='sample1 sample2',
#                       help='User defined labels instead of default labels from '
#                       'file names. '
#                       'Multiple labels have to be separated by spaces, e.g., '
#                       '--labels sample1 sample2 sample3',
#                       nargs='+')

args = parser.parse_args()

# #from pybedtools.helpers import get_chromsizes_from_ucsc
# #chromsizes = get_chromsizes_from_ucsc('hg19')
# intervals = BedTool().window_maker(w=10000, genome='hg19')
# BedTool().annotate(
#     counts=True,
#     i=intervals,
#     files=['events/uvCLAP_run1_hg19_F1_hDHX9_repA.bed.gz', 'events/uvCLAP_run1_hg19_F1_hDHX9_repB.bed.gz'])
