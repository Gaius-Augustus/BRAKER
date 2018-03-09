wd=test4

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=genome.fa --prot_seq=prot.fa --prg=gth --trainFromGth --workingdir=$wd ) &> test4.log
