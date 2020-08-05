wd=test8_restart3
oldDir=test8

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test8.sh before running test8_restart3.sh!"
else
    species=$(cat $oldDir/braker.log | perl -ne 'if(m/AUGUSTUS parameter set with name ([^.]+)\./){print $1;}')
    ( time braker.pl --genome=../genome.fa --esmode --skipAllTraining --softmasking --workingdir=$wd --species=$species --cores 8 ) &> test8_restart3.log
fi
