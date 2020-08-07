wd=test3_restart1
oldDir=test3

if [ -d $wd ]; then
    rm -r $wd
fi

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test3.sh before running test3_restart1.sh!"
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --etpmode --softmasking --workingdir=$wd --cores 8 --gm_max_intergenic 10000) &> test3_restart1.log
fi
