wd=test8

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=genome.fa --esmode --softmasking --workingdir=$wd --cores=2 ) &> test8.log
