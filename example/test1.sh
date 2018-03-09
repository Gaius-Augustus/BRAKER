wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=genome.fa --bam=RNAseq.bam --workingdir=$wd ) &> test1.log
