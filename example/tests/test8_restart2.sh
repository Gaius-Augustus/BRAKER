wd=test8_restart2
oldDir=test8

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test8.sh before running test8_restart2.sh!"
else
    ( time braker.pl --genome=../genome.fa --esmode --geneMarkGtf=$oldDir/GeneMark-ES/genemark.gtf --softmasking --workingdir=$wd --cores 8 ) &> test8_restart2.log
fi
