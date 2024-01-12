wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# Note:
# The file ../RNASeq.bam is not contained in the github repository!
# Make sure that you downloaded this file with
#    wget http://topaz.gatech.edu/GeneMark/Braker/RNAseq.bam
# before running this test.
# Alternatively, you can replace --bam=../RNASeq.bam by --hints=../RNASeq.hints

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes. Also
# --skipOptimize is used here only to make the test run faster.
# Use an appropriate BUSCO lineage in real life use case. eukaryota_odb10 is used here
# only to make the test run faster.

export GENEMARK_PATH=$GENEMARK_PATH/gmes

( time braker.pl --genome=../genome.fa --bam=../RNAseq.bam --workingdir=$wd --threads 11 --gm_max_intergenic 10000 --skipOptimize --busco_lineage eukaryota_odb10 ) &> test1.log
