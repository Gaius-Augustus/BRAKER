wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome.fa --hints=../prothint_augustus.gff,../evidence_augustus.gff --prothints=../prothint.gff --evidence=../evidence.gff --bam=../RNAseq.bam --etpmode --softmasking --workingdir=$wd ) &> test3.log
