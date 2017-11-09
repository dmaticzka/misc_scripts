#!/urs/bin/env bash
[ $# -ne 4 ] &&  echo "Usage:  $(basename $0) pos_dir neg_dir work_dir genometab" && echo "Instead of: $(basename $0) $*"  && exit 1

JAMM=JAMM.sh
JAMMTHREADS=23

POSDIR=$1
NEGDIR=$2
WORKDIR=$3
GENOMETAB=$4

rm -rf $WORKDIR

# split up positives into plus and minus strand
# and convert to BEDPE
mkdir -p $WORKDIR/pos/plus $WORKDIR/pos/minus
for BED in $POSDIR/*.bed
do
    cat $BED | \
    bedtools sort | \
	awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1, $2, $3, $4, $5, $6, $6}' | \
    awk '$9=="+"' > $WORKDIR/pos/plus/`basename $BED`
    cat $BED | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1, $2, $3, $4, $5, $6, $6}' | \
    awk '$9=="-"' > $WORKDIR/pos/minus/`basename $BED`
done

# split up negatives into plus and minus
mkdir -p $WORKDIR/neg/plus $WORKDIR/neg/minus
for BED in $NEGDIR/*.bed
do
    cat $BED | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1, $2, $3, $4, $5, $6, $6}' | \
    awk '$9=="+"' > $WORKDIR/neg/plus/`basename $BED`
    cat $BED | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $1, $2, $3, $4, $5, $6, $6}' | \
    awk '$9=="-"' > $WORKDIR/neg/minus/`basename $BED`
done

# run jamm
## FIXBELOW
for STRAND in plus minus
do
    mkdir -p $WORKDIR/$STRAND
    time $JAMM \
    -s $WORKDIR/pos/$STRAND \
    -c $WORKDIR/neg/$STRAND \
    -g $GENOMETAB \
    -o $WORKDIR/$STRAND \
    -d y \
    -t paired \
    -b 50 \
    -w 1 \
    -m normal \
    -p $JAMMTHREADS 2>&1 | tee $WORKDIR/jamm_$STRAND.log; \
done

# combine strands
# collect data
( cat $WORKDIR/plus/peaks/all.peaks.narrowPeak | awk '{OFS="\t"}{$6="+"; print $0}'; \
 cat $WORKDIR/minus/peaks/all.peaks.narrowPeak | awk '{OFS="\t"}{$6="-"; print $0}'; ) | \
 bedtools sort > $WORKDIR/all.peaks.narrowPeak;
( cat $WORKDIR/plus/peaks/filtered.peaks.narrowPeak | awk '{OFS="\t"}{$6="+"; print $0}'; \
 cat $WORKDIR/minus/peaks/filtered.peaks.narrowPeak | awk '{OFS="\t"}{$6="-"; print $0}'; ) | \
 bedtools sort > $WORKDIR/filtered.peaks.narrowPeak;
