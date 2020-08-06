wd=test8

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

( time braker.pl --genome=../genome.fa --esmode --softmasking --workingdir=$wd --cores 8 --gm_max_intergenic 10000 ) &> test8.log
