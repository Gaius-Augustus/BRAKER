wd=test2

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=../genome.fa --prot_seq=../orthodb_small.fa --epmode --softmasking --workingdir=$wd --cores=11 ) &> test2.log
