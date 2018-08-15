wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=genome.fa --hints=ep.hints --bam=RNAseq.bam --etpmode --softmasking --workingdir=$wd ) &> test3.log
