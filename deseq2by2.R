#! /usr/bin/env Rscript

require(argparse)
require(plyr)
require(dplyr)
require(tidyr)

parser <- ArgumentParser(
    description='run deseq2 for two positives and two controls that allows to manually set library normalization factors')
parser$add_argument(
    'pos',
    metavar='pos_counts',
    type='character',
    nargs=2,
    help="two tables containing counts of positive replicates")
parser$add_argument(
    'neg',
    metavar='neg_counts',
    type='character',
    nargs=2,
    help="two tables containing counts of control replicates")
parser$add_argument(
    '-n',
    '--normalization-factors',
    type='double',
    nargs='+',
    help="if given, these values are used for library normalization")

args <- parser$parse_args()

args$normalization_factors
if (!is.null(args$normalization_factors)) {
    if (length(args$normalization_factors) != 4) {
        stop('must supply exactly four factors when doing manual library normalization')
    }
}
