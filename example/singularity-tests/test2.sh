#!/bin/bash

# Author: Katharina J. hoff

# Contact: katharina.hoff@uni-greifswald.de

# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
# Then run "bash test2.sh".

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh ."
    echo "Then execute on the host with \"bash test2.sh\"".
    exit 1
fi

# remove output directory if it already exists

wd=test2

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=/opt/BRAKER/example/genome.fa --prot_seq=/opt/BRAKER/example/proteins.fa --workingdir=${wd} --threads 8 --gm_max_intergenic 10000 --skipOptimize --busco_lineage eukaryota_odb10 &> test2.log
            # Important: the options --gm_max_intergenic 10000 --skipOptimize should never be applied to a real life run!!!
            # They were only introduced to speed up the test. Please delete them from the script if you use it for real data analysis.


