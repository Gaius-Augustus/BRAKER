wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

 time braker.pl --genome=../genome.fa --bam=../RNAseq.bam --softmasking --workingdir=$wd --cleanup ) &> test1.log
