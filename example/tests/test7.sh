wd=test7

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=../genome.fa --bam=../RNAseq.bam --species=fly --skipAllTraining --softmasking --workingdir=$wd --cleanup ) &> test7.log
