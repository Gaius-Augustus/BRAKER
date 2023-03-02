wd=test2_restart3
oldDir=test2

if [ -d $wd ]; then
    rm -r $wd
fi

export GENEMARK_PATH=$GENEMARK_PATH/gmes

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test2.sh before running test2_restart2.sh!"
else
    species=$(cat $oldDir/braker.log | perl -ne 'if(m/AUGUSTUS parameter set with name ([^.]+)\./){print $1;}')
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --skipAllTraining --species=$species --workingdir=$wd --threads 8 --skipOptimize ) &> test2_restart3.log
fi
