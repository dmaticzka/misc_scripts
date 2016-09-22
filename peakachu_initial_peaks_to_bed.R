#! /usr/bin/env Rscript

suppressPackageStartupMessages({
    library('argparse')
    library('dplyr')
})

parser <- ArgumentParser(
    description='convert PEAKachu csv output of potential peaks to bed (PEAKachu adaptive using DESeq2).')
parser$add_argument(
    'file',
    type='character',
    nargs=1,
    help='PEAKachu adaptive DESeq2 output (initial_peaks.csv).')
args <- parser$parse_args()

outfile <- paste0(args$file,'.bed')

read.table(args$file, head=T) %>%
    select(replicon, peak_start, peak_end, log2FoldChange, padj, peak_strand) %>%
    mutate(log2FoldChange=paste0('l2fc', log2FoldChange, ';padj:', padj), padj=255) %>%
    write.table(outfile, sep='\t', quote=F, row.names=F, col.names=F)
