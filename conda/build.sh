mkdir -p $PREFIX/bin

#PROGS="multiBedSummary.py peakachu_initial_peaks_to_bed.R prepare_graphprot_seqs.py runjamm_from_bed6.sh simple_deseq2.R yass_reverse_repeats_to_bed12.R"

# manually distribute files
#for PROG in $PROGS
#do
#    cp -v $SRC_DIR/$PROG $PREFIX/bin/
#done

# copy all
cp -v $SRC_DIR/* $PREFIX/bin/
