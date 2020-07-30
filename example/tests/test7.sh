wd=test7

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=../genome2.fa --bam=../RNAseq2.bam --species=arabidopsis --skipAllTraining --softmasking --workingdir=$wd --cores 8 ) &> test7.log
