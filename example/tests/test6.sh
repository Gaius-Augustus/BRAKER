wd=test6

if [ -d $wd ]; then
    rm -r $wd
fi

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

( time braker.pl --genome=../genome2.fa --prot_seq=../proteins2.fa --prg=gth --bam=../RNAseq2.bam --gth2traingenes --softmasking --workingdir=$wd --cores 8 --gm_max_intergenic 10000 ) &> test6.log
