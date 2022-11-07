wd=test5

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

export GENEMARK_PATH=$GENEMARK_PATH/gmes

( time braker.pl --genome=../genome.fa --esmode --workingdir=$wd --threads 8 --gm_max_intergenic 10000 --skipOptimize ) &> test5.log
