#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))

# parse args
parser <- ArgumentParser(
    description='run deseq2 for two positives and two controls. allows to manually set library normalization factors. count tables are expected to contain a header line as e.g. supplied by featureCounts and have ids in the first and counts in the last column.')
parser$add_argument(
    'pos',
    type='character',
    nargs=2,
    help='two tables containing counts of positive replicates.')
parser$add_argument(
    'neg',
    type='character',
    nargs=2,
    help='two tables containing counts of control replicates.')
parser$add_argument(
    '-n',
    '--normalization-factors',
    type='double',
    nargs='+',
    help='if given, these values are used for library normalization. expects 4 space separated values in order pos1, pos2, neg1, neg2.')
parser$add_argument(
    '-p',
    '--prefix',
    type='character',
    default='results_deseq2by2',
    help='prefix output files with this string (default: results_deseq2by2)')
args <- parser$parse_args()

# check normalization factors
if (!is.null(args$normalization_factors)) {
    if (length(args$normalization_factors) != 4) {
        stop('must supply exactly four factors when doing manual library normalization')
    }
}

# load counts
pos1 <- read.table(args$pos[1], head=T) %>% select(1,ncol(.))
pos2 <- read.table(args$pos[2], head=T) %>% select(1,ncol(.))
neg1 <- read.table(args$neg[1], head=T) %>% select(1,ncol(.))
neg2 <- read.table(args$neg[2], head=T) %>% select(1,ncol(.))

# prepare counts for deseq2
d <- join_all(list(pos1,pos2,neg1,neg2))
rownames(d) <- d[,1]
d[,1] <- NULL

# prepare sample table
sampleTable <-
    data.frame(
        condition=c('pos','pos','neg','neg'),
        colnames(d))

print('using sample table:')
sampleTable
print('')
print('')

# create DESeq2 object
dds <-
    DESeqDataSetFromMatrix(
        countData = d,
        colData = sampleTable,
        design = ~ condition)
# pre-filtering to reduce memory size
dds <- dds[ rowSums(counts(dds)) > 1, ]

# library normalization
dds <- estimateSizeFactors(dds)
print('size factors estimated by deseq2:')
print(dds$sizeFactor)
if (!is.null(args$normalization_factors)) {
    dds$sizeFactor <- args$normalization_factors
    print('manually set size factors:')
    print(dds$sizeFactor)
}

# finish deseq2 analysis
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

# some output
summary(res)

# convert to data frame with separate id column
dres <- data.frame(res) %>% add_rownames(var='id')

# write tables full and significant to files
resFull <- arrange(dres, padj)
write.table(
    resFull,
    paste0(args$prefix,'.deseq2_full.csv'),
    sep="\t",
    quote=F,
    row.names=F)
resSig <- filter(resFull, padj < 0.1)
write.table(
    resSig,
    sep="\t",
    paste0(args$prefix,'.deseq2_significant.csv'),
    quote=F,
    row.names=F)
