# Copy this script into the folder where you want to execute it prior modification and
# prior running it, e.g.:
# singularity exec -B $PWD:$PWD /opt/BRAKER/example/singularity-tests/test1.sh .

# Next, modifiy the ETP folder time to your system
ETP=${HOME}/ETP

wd=test1

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes.

singularity exec -B $PWD:$PWD braker3.sif braker.pl --genome=/opt/BRAKER/example/genome.fa --bam=/opt/BRAKER/example/RNAseq.bam --softmasking --workingdir=${wd} --threads 8 --GENEMARK_PATH=${ETP}/bin/gmes \
	    --gm_max_intergenic 10000 --skipOptimize


