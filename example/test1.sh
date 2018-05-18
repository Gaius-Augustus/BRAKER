wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=genome.fa --bam=RNAseq.bam --workingdir=$wd --UTR=on --cores=6 --AUGUSTUS_ab_initio --softmasking --skipOptimize ) &> test1.log
