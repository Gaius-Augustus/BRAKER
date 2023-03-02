wd=test3_restart2
oldDir=test3

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test3.sh before running test3_restart2.sh!"
else
    ( time braker.pl --genome=../genome.fa --gmetp_results_dir=$oldDir/GeneMark-ETP/ --workingdir=$wd --threads 8 --skipOptimize ) &> test3_restart2.log
fi
