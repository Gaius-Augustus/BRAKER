wd=test10

if [ -d $wd ]; then
    rm -r $wd
fi

# if you have not retrieved the data, uncomment the following lines for a one time execution:
# cd ..
# wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
# tar xvf odb10_arthropoda_fasta.tar.gz
# cat arthropoda/Rawdata/* > proteins.fasta
# sed -i "s/\.//" proteins.fasta
# rm -r  arthropoda
# cd tests

( time braker.pl --genome=../genome.fa --prot_seq=../proteins.fasta --prg=ph --epmode --softmasking --workingdir=$wd --ALIGNMENT_TOOL_PATH=/home/katharina/git/ProtHint/bin ) &> test10.log
