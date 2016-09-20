#! /usr/bin/env Rscript

suppressPackageStartupMessages({
    library('argparse')
    library('tidyr')
    library('plyr')
    library('dplyr')
    library('DESeq2')
    library('pheatmap')
    library('RColorBrewer')
    library('ggplot2')
})

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
    help='prefix output files with this string (default: results_deseq2by2).')
parser$add_argument(
    '-s',
    '--significance',
    type='double',
    default='0.1',
    help='set threshold for adjusted p-values written to prefix.deseq2_significant.csv (default: 0.1).')
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
        source=c(args$pos, args$neg),
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
res <- results(dds, contrast=c('condition','pos','neg'))

# some output
summary(res)

# convert to data frame with separate id column
dres <- data.frame(res) %>% add_rownames(var='id')

# write tables full and significant to files
resFull <- arrange(dres, padj)
write.table(
    resFull,
    paste0(args$prefix,'.deseq2_full.csv'),
    sep='\t',
    quote=F,
    row.names=F)
resSig <- filter(resFull, padj < args$significance)
write.table(
    resSig,
    sep='\t',
    paste0(args$prefix,'.deseq2_significant.csv'),
    quote=F,
    row.names=F)

# plot clustering via sample-to-sample distances
plot_sample_dists <- function(values) {
    sampleDists <- dist(t(assay(values)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(values$condition, values$source, sep='-')
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
    pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors)
}

# do plots
pdf(paste0(args$prefix,'.deseq2_plots.pdf'), paper='a4')

vts <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
plot_sample_dists(dds)
plot_sample_dists(rld)
plot_sample_dists(vts)

data <- plotPCA(rld, intgroup=c('condition', 'source'), returnData=TRUE)
percentVar <- round(100 * attr(data, 'percentVar'))
ggplot(data, aes(PC1, PC2, color=condition, shape=source)) +
geom_point(size=3) +
xlab(paste0('PC1: ',percentVar[1],'% variance')) +
ylab(paste0('PC2: ',percentVar[2],'% variance'))

# # need meaningful labels, because from Galaxy, sample names are random
# labs <- paste0(seq_len(ncol(dds)), ': ', do.call(paste, as.list(colData(dds)[factors])))
# dat <- assay(rld)
# colnames(dat) <- labs
# distsRL <- dist(t(dat))
# mat <- as.matrix(distsRL)
# hc <- hclust(distsRL)
# hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
# heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace='none', col = rev(hmcol),
#           main='Sample-to-sample distances', margin=c(13,13))

plotDispEsts(dds, main='Dispersion estimates')

dev.off()
