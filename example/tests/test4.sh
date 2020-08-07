wd=test4

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~7 minutes. The fast runtime of
# this test is mostly caused by generating a low number of training genes.
# Note that this approach does not scale well with increasing genome size
# and the number of proteins in a protein database. The runtime on a full
# genome will be much slower than with the command used in test2.sh.


( time braker.pl --genome=../genome.fa --prot_seq=../proteins.fa --prg=gth --trainFromGth --softmasking --workingdir=$wd --cores 8) &> test4.log
