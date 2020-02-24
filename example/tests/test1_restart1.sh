wd=test1_restart1
oldDir=test1

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test1.sh before running test1_restart1.sh!"  
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --softmasking --workingdir=$wd ) &> test1_restart1.log
fi
