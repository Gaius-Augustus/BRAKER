wd=test1_restart2
oldDir=test1

if [ -d $wd ]; then
    rm -r $wd
fi

export GENEMARK_PATH=$GENEMARK_PATH/gmes

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test1.sh before running test1_restart2.sh!"
else
    ( time braker.pl --genome=../genome.fa --hints=$oldDir/hintsfile.gff --geneMarkGtf=$oldDir/GeneMark-ET/genemark.gtf --workingdir=$wd --threads 8 --skipOptimize ) &> test1_restart2.log
fi
