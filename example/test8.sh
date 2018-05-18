wd=test8

if [ -d $wd ]; then
    rm -r $wd
fi

( time braker.pl --genome=genome.fa --workingdir=$wd --cores=6 --softmasking --skipOptimize --esmode --cores=8 --skipOptimize ) &> test8.log
