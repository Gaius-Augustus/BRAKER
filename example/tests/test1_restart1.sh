wd=test1_restart1
oldDir=test1

if [ -d $wd ]; then
    rm -r $wd
fi

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

export GENEMARK_PATH=$GENEMARK_PATH/gmes

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test1.sh before running test1_restart1.sh!"
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --workingdir=$wd --threads 8 --gm_max_intergenic 10000 --skipOptimize ) &> test1_restart1.log
fi
