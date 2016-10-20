#! /usr/bin/env Rscript

suppressPackageStartupMessages({
    library('argparse')
    library('tidyr')
    library('dplyr')})

parser <- ArgumentParser(
    description='load yass output in bed format (sans track lines) and convert to bed12 format for display where repeat regions are the two blocks. fasta ids should contain the genomic coordinates of the sequence in the format as given by bedtools getfasta (eg "chr1:12057296-12061695(+)").')

parser$add_argument(
    '-i',
    '--yassbed',
    type='character',
    help='yass output in bed format')

parser$add_argument(
    '-o',
    '--outbed',
    type='character',
    help='bed12 output')

args <- parser$parse_args()

print(args$yassbed)
print(args$outbed)

# read yass bed file
yb <- read.table(
    args$yassbed,
    col.names=c(
        "bedid",
        "s_start",
        "s_end",
        "strand",
        "q_start",
        "q_end",
        "bitscore",
        "evalue",
        "identity"))

# prepare data for processing
yb <- yb %>%
    # split first column containing absolute coordinates
    separate(
        bedid,
        c("chr","chr_range","chr_strand","nothing"),
        sep = "[^[:alnum:]^+^-]+",
        convert=T,
        remove=F) %>%
    # second split to get chromosomal start stop coordinates separated by -
    separate(
        chr_range,
        c("chr_start","chr_end"),
        convert=T)

# generate bed 12 for display
d_bed12 <- yb %>%
    # only keep s regions >= 150
    filter(s_end-s_start>=150) %>%
    rowwise %>%
    # get outer coordinates
    mutate(
        start=chr_start+q_start,
        end=chr_start+s_end) %>%
    # get block sizes
    mutate(
        bsize_1=q_end-q_start,
        bsize_2=s_end-s_start)

nreverse <- d_bed12 %>% filter(s_start<q_start) %>% nrow
if (nreverse != 0) {
    stop("Error, encountered pair that has s left of q.")
}

# get actual bed12 format
d_bed12 <- d_bed12 %>%
    transmute(
        chr=chr,
        start=start,
        end=end,
        name=paste0(
            "qlength:",
            q_end-q_start,
            ",distance:",
            s_start-q_end,
            ",identity:",
            identity),
        score=q_end-q_start,
        strand=chr_strand,
        thickStart=start,
        thickEnd=end,
        itemRGB="0,0,255",
        blockCount=2,
        blockSizes=paste(bsize_1,bsize_2,sep=","),
        blockStarts=paste(0,s_start-q_start,sep=","))

# write bed12 to disk
write.table(
    d_bed12,
    args$outbed,
    sep="\t",
    row.names=F,
    col.names=F,
    quote=F)
