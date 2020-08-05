wd=test4

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome.fa --prot_seq=../proteins.fa --prg=gth --trainFromGth --softmasking --workingdir=$wd --cores 8) &> test4.log
