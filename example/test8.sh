wd=test8

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=genome.fa --esmode --softmasking --workingdir=$wd ) &> test8.log
