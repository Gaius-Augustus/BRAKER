wd=test7_restart1
oldDir=test7

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test7.sh before running test7_restart1.sh!"
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --species=arabidopsis --skipAllTraining --softmasking --workingdir=$wd --cores 8 ) &> test7_restart1.log
fi

