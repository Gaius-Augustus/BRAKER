wd=test6

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome.fa --prot_seq=../prot.fa --prg=gth --bam=../RNAseq.bam --gth2traingenes --softmasking --workingdir=$wd --cleanup ) &> test6.log
