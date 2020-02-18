wd=test2_restart1
oldDir=test2

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test2.sh before running test2_restart1.sh!"  
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --evidence=$oldDir/evidence.gff --prothints=$oldDir/prothint.gff --epmode --softmasking --workingdir=$wd --cores=11 ) &> test2_restart1.log
fi
