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
    library('gplots')
})

# parse args
parser <- ArgumentParser(
    description='run deseq2 with a simple positive/negative contrast. allows to manually set library normalization factors. count tables are expected to contain a header line as e.g. supplied by featureCounts and have ids in the first and counts in the last column. some summary plots are written to pdf.')
parser$add_argument(
    '--pos',
    type='character',
    nargs='+',
    help='two tables containing counts of positive replicates.')
parser$add_argument(
    '--neg',
    type='character',
    nargs='+',
    help='two tables containing counts of control replicates.')
parser$add_argument(
    '-n',
    '--normalization-factors',
    type='double',
    nargs='+',
    help='if given, these values are used for library normalization. expects space separated values in order pos_1..pos_n, neg_1, neg_n.')
parser$add_argument(
    '-p',
    '--prefix',
    type='character',
    default='simple_deseq2',
    help='prefix output files with this string (default: simple_deseq2).')
parser$add_argument(
    '-s',
    '--significance',
    type='double',
    default='0.1',
    help='set threshold for adjusted p-values written to prefix.deseq2_significant.csv (default: 0.1).')
args <- parser$parse_args()

# check normalization factors
if (!is.null(args$normalization_factors)) {
    if (length(args$normalization_factors) != length(args$pos) + length(args$neg)) {
        stop('must supply as many normalization factors as replicates library normalization')
    }
}

# load counts
load_counts <- function(fn) {
    return(read.table(fn, head=T) %>% select(1,ncol(.)))
}
countfiles <- c(args$pos, args$neg)
d <- lapply(countfiles, load_counts)
names(d) <- countfiles

# prepare counts for deseq2
d <- join_all(d)
rownames(d) <- d[,1]
d[,1] <- NULL

# prepare sample table
sampleTable <-
    data.frame(
        source=c(args$pos, args$neg),
        condition=c(
            rep('pos', length(args$pos)),
            rep('neg', length(args$neg))),
        colnames(d))

print('using sample table:')
sampleTable
print('')
print('')

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
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
plot_sample_dists <- function(values, title) {
    sampleDists <- dist(t(assay(values)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(values$condition, values$source, sep='-')
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
    pheatmap(
        sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors,
        main=title)
}

# do plots
pdf(paste0(args$prefix,'.deseq2_plots.pdf'), paper='a4')

# MA-plot
plotMA(res, main="DESeq2")

vts <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# plot sample-to sample distances
plot_sample_dists(dds, title="Sample-to-sample distances (untransformed)")
plot_sample_dists(rld, title="Sample-to-sample distances (regularized log)")
plot_sample_dists(vts, title="Sample-to-sample distances (variance stabilizing transformation)")

# plot PCA
data <- plotPCA(rld, intgroup=c('condition', 'source'), returnData=TRUE)
percentVar <- round(100 * attr(data, 'percentVar'))
ggplot(data, aes(PC1, PC2, color=condition, shape=source)) +
geom_point(size=3) +
xlab(paste0('PC1: ',percentVar[1],'% variance')) +
ylab(paste0('PC2: ',percentVar[2],'% variance'))

# plot dispersion estimates
plotDispEsts(dds, main='Dispersion estimates')

dev.off()
