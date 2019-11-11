wd=test2

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=../genome.fa --hints=../prothint_augustus.gff,../evidence_augustus.gff --evidence=../evidence.gff --prothints=../prothint.gff --epmode --softmasking --workingdir=$wd ) &> test2.log
