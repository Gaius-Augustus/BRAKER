wd=test3_4

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# Note:
# The file ../RNASeq.bam is not contained in the github repository!
# Make sure that you downloaded this file with
#     wget http://topaz.gatech.edu/GeneMark/Braker/RNAseq.bam
# before running this test.
# Alternatively, you can replace --bam=../RNASeq.bam by --hints=../RNASeq.hints

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

# For instructions on how to prepare the proteins.fa input file from OrthoDB,
# see https://github.com/gatech-genemark/ProtHint#protein-database-preparation

( time braker.pl --genome=../genome.fa --prot_seq=../proteins.fa --rnaseq_sets_ids=ERR5767212 --workingdir=$wd --threads 8 ) &> test3_4.log
