wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome.fa --prot_seq=../orthodb_small.fa --bam=../RNAseq.bam --etpmode --softmasking --workingdir=$wd ) &> test3.log
