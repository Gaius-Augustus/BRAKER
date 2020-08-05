wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi

# Note:
# The file ../RNASeq2.bam is not contained in the github repository!
# Make sure that you downloaded this file with
# wget http://topaz.gatech.edu/GeneMark/Braker/RNAseq2.bam
# before running this test.
# Alternatively, you can replace --bam=../RNASeq.bam by --hints=../RNASeq.hints

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

( time braker.pl --genome=../genome2.fa --prot_seq=../proteins2.fa --bam=../RNAseq2.bam --etpmode --softmasking --workingdir=$wd --cores 8 --gm_max_intergenic 10000 ) &> test3.log
