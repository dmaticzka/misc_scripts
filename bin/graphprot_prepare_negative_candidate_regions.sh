#!/usr/bin/env bash

[ $# -ne 4 ] &&
echo "Usage:  $(basename $0) PEAKS GENEDEF GENOME_LIMITS OUTPUT instead of: $(basename $0) $*"  && exit 1

PEAKS=$1
GENEDEF=$2
GENOME_LIMITS=$3
NEGATIVE_CANDIDATE_REGIONS=$4

# select targeted genes
bedtools sort -i $PEAKS |
bedtools intersect -s -u \
-a $GENEDEF \
-b $PEAKS |
cut -f 1-6 > $PEAKS.targets

# remove peaks and surrounding regions from targeted genes
bedtools sort -i $PEAKS |
bedtools slop -b 100 -i - -g $GENOME_LIMITS |
bedtools subtract -s -a $PEAKS.targets -b - > $NEGATIVE_CANDIDATE_REGIONS
