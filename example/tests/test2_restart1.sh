wd=test2_restart1
oldDir=test2

if [ -d $wd ]; then
    rm -r $wd
fi

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

# Warning: the result of this particular test can slightly differ from the full
# test run (test2.sh) because the protein hints supplied to GeneMark-EP+ in this
# test are from ProtHint's second iteration, while in the regular run, GeneMark-EP+
# uses result from ProtHint's first iteration (see BRAKER2 paper for details).

export GENEMARK_PATH=$GENEMARK_PATH/gmes

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test2.sh before running test2_restart1.sh!"
else
  ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --workingdir=$wd --threads 8 --gm_max_intergenic 10000 --skipOptimize ) &> test2_restart1.log
fi
