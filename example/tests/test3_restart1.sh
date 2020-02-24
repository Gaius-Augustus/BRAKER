wd=test3_restart1
oldDir=test3

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test3.sh before running test3_restart1.sh!"  
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --evidence=$oldDir/evidence.gff --genemark_hintsfile=$oldDir/genemark_hintsfile.gff --etpmode --softmasking --workingdir=$wd ) &> test3_restart1.log
fi
