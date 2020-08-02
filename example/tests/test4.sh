wd=test4

if [ -d $wd ]; then
    rm -r $wd
fi


( time braker.pl --genome=../genome2.fa --prot_seq=../proteins2.fa --prg=gth --trainFromGth --softmasking --workingdir=$wd --cores 8) &> test4.log
