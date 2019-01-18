wd=test5

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome.fa --prot_seq=../prot.fa --prg=gth --bam=../RNAseq.bam --softmasking --workingdir=$wd --cleanup ) &> test5.log
