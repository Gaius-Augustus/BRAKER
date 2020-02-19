wd=test5_restart1

oldDir=test5

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test5.sh before running test5_restart1.sh!"  
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --softmasking --workingdir=$wd --cores=11 ) &> test5_restart1.log
fi


