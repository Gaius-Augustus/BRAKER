wd=test4

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~1 minute.

# Note:
# The file ../RNASeq.bam is not contained in the github repository!
# Make sure that you downloaded this file with
#    wget http://topaz.gatech.edu/GeneMark/Braker/RNAseq.bam
# before running this test.
# Alternatively, you can replace --bam=../RNASeq.bam by --hints=../RNASeq.hints

( time braker.pl --genome=../genome.fa --bam=../RNAseq.bam --species=arabidopsis --skipAllTraining --workingdir=$wd --threads 8 ) &> test4.log
