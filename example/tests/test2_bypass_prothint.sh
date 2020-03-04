wd=test2_bypass_prothint
oldDir=test2

if [ -d $wd ]; then
    rm -r $wd
fi

if [ ! -d $oldDir ] ; then
  echo "ERROR: Directory (with contents) of old BRAKER run $oldDir does not exist, yet. Please run test2.sh after adding --nocleanup, and before running test2_bypass_prothint.sh!"  
else
    # files will be written to BRAKER/example folder
    log_reg_prothints.pl --prothint=$oldDir/prothint.gff --out=../prothint_reg_augustus.gff
    cat $oldDir/evidence_augustus.gff | perl -pe 's/src=P/src=M/' > ../evidence_manual.gff
    cat $oldDir/top_chains.gff | perl -pe 's/src=P/src=C/;' > ../top_chains_c.gff
    ( time braker.pl --genome=../genome.fa --hints=../prothint_reg_augustus.gff,../evidence_manual.gff,../top_chains_c.gff --evidence=$oldDir/evidence.gff --prothints=$oldDir/prothint.gff --epmode --softmasking --workingdir=$wd ) &> test2_bypass_prothint.log
fi
