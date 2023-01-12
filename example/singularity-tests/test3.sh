#!/bin/bash

# Author: Katharina J. hoff
# Contact: katharina.hoff@uni-greifswald.de
# Date: Jan 12th 2023

# Copy this script into the folder where you want to execute it, e.g.:
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# Then run "bash test3.sh".

# Check wheter ETP is installed

if [[ -z "${ETP}" ]]; then
    echo ""
    echo "Variable ETP is undefined."
    echo "You need to install GeneMark-ETP in the host home directory before running the BRAKER container."
    echo "After installation (get instructions with \"singularity exec braker3.sif print_braker3_setup.py\")"
    echo "export the ETP environment variable as follows:"
    echo ""
    echo "export ETP=\${HOME}/ETP_folder/bin"
    echo ""
    echo "(You will have to modify the folder name to where the etp_release.pl script sits,"
    echo "that might not be \${HOME}/ETP_folder on your system.)"
    echo ""
    exit 1      
else
    if [ ! -f "${ETP}/etp_release.pl" ]; then
	echo "Could not find ${ETP}/etp_release.pl. GeneMark-ETP has not been installed, correctly."
	echo "Get help with \"singularity exec braker3.sif print_braker3_setup.p\""
	exit 1
    fi
    if [ ! -f "${HOME}/.gm_key" ]; then
	echo "GeneMark-ETP key not found at ${HOME}/.gm_key ."
	echo "It is not impossible that things are ok. However, most likely, something is wrong."
	echo "Get helo with \"singularity exec braker3.sif print_braker3_setup.p\" if needed."
    fi
fi

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
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh ."
    echo "Then execute on the host with \"bash test3.sh\"".
    exit 1
fi

# remove output directory if it already exists

wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=/opt/BRAKER/example/genome.fa --prot_seq=/opt/BRAKER/example/proteins.fa --bam=/opt/BRAKER/example/RNAseq.bam --softmasking --workingdir=${wd} \
	    --GENEMARK_PATH=${ETP} --PROTHINT_PATH=${ETP}/gmes/ProtHint/bin --threads 8 --gm_max_intergenic 10000 --skipOptimize

	    # Important: the options --gm_max_intergenic 10000 --skipOptimize should never be applied to a real life run!!!                                   
            # They were only introduced to speed up the test. Please delete them from the script if you use it for real data analysis. 
