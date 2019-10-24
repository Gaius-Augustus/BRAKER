wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

# Note:
# The file ../RNASeq.bam is not contained in the github repository!
# Make sure that you downloaded this file with
# wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam
# before running this test.
# Alternatively, you can replace --bam=../RNASeq.bam by --hints=../RNASeq.hints

( time braker.pl --genome=../genome.fa --bam=../RNAseq.bam --softmasking --workingdir=$wd ) &> test1.log
