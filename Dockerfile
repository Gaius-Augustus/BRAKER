# Distributed under the terms of the Modified BSD License.
ARG OWNER=jupyter
ARG BASE_CONTAINER=$OWNER/minimal-notebook
FROM $BASE_CONTAINER as base

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    # for cython: https://cython.readthedocs.io/en/latest/src/quickstart/install.html
    build-essential \
    # for latex labels
    cm-super \
    dvipng \
    # for matplotlib anim
    ffmpeg \
    time && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN  apt update && \
     apt-get install -y --no-install-recommends \
     man-db \
     g++ \
     less \
     zlib1g-dev \
     && \
     apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cd /opt && \ 
    git clone --recursive https://github.com/clwgg/seqstats && \
    cd seqstats && \
    make 
    
ENV PATH=${PATH}:/opt/seqstats

# cdbfasta
RUN cd /opt && \ 
    git clone https://github.com/gpertea/cdbfasta.git && \
    cd cdbfasta && \
    make 
    
ENV PATH=${PATH}:/opt/cdbfasta

# tsebra
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/TSEBRA

ENV PATH=${PATH}:/opt/TSEBRA/bin

# makehub
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/MakeHub.git && \
    cd MakeHub && \
    git checkout braker3 && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredCheck && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/faToTwoBit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/gtfToGenePred && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/hgGcPercent && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/ixIxx && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitInfo && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBed && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/genePredToBigGenePred && \
    chmod u+x bedToBigBed genePredCheck faToTwoBit gtfToGenePred hgGcPercent ixIxx  twoBitInfo wigToBigWig genePredToBed genePredToBigGenePred make_hub.py

ENV PATH=${PATH}:/opt/MakeHub

# augustus

ENV AUGUSTUS_CONFIG_PATH=/usr/share/augustus/config/
ENV AUGUSTUS_BIN_PATH=/usr/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts/

RUN apt update && \ 
    apt install -yq  augustus augustus-data augustus-doc && \
    apt clean all && \
    fix-permissions "${AUGUSTUS_CONFIG_PATH}"

# perl dependencies of BRAKER and GeneMark-ETP+

RUn apt update && \
    apt install -yq libyaml-perl libhash-merge-perl libparallel-forkmanager-perl libscalar-util-numeric-perl libclass-data-inheritable-perl libexception-class-perl libtest-pod-perl libfile-which-perl libmce-perl libthread-queue-perl libmath-utils-perl libscalar-list-utils-perl && \
    apt clean all
    
# patch Augustus scripts (because Debian package is often outdated, this way we never need to worry)
RUN cd /usr/share/augustus/scripts && \
    rm optimize_augustus.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/optimize_augustus.pl && \
    rm aa2nonred.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/aa2nonred.pl && \
    rm gff2gbSmallDNA.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/gff2gbSmallDNA.pl && \
    rm new_species.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/new_species.pl && \
    rm filterGenesIn.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/filterGenesIn.pl && \
    rm filterGenesIn_mRNAname.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/filterGenesIn_mRNAname.pl && \
    rm filterGenes.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/filterGenes.pl && \
    rm join_mult_hints.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/join_mult_hints.pl && \
    rm randomSplit.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/randomSplit.pl && \
    rm join_aug_pred.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/join_aug_pred.pl && \
    rm getAnnoFastaFromJoingenes.py && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/getAnnoFastaFromJoingenes.py && \
    rm gtf2gff.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/gtf2gff.pl && \
    rm splitMfasta.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/splitMfasta.pl && \
    rm createAugustusJoblist.pl && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/scripts/createAugustusJoblist.pl && \
    chmod a+x optimize_augustus.pl aa2nonred.pl gff2gbSmallDNA.pl new_species.pl filterGenesIn_mRNAname.pl filterGenes.pl filterGenesIn.pl join_mult_hints.pl randomSplit.pl join_aug_pred.pl getAnnoFastaFromJoingenes.py gtf2gff.pl splitMfasta.pl createAugustusJoblist.pl

RUN cd /usr/share/augustus/config && \
    mkdir parameters && \
    cd parameters && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/config/parameters/AUG_CMDLN_PARAMETERS.md && \
    wget https://raw.githubusercontent.com/Gaius-Augustus/Augustus/master/config/parameters/aug_cmdln_parameters.json && \
    chmod a+r AUG_CMDLN_PARAMETERS.md aug_cmdln_parameters.json && \
    cd .. && \
    chmod a+r parameters && \
    chmod a+x parameters

# bamtools (ETP+)

RUN apt update && \
    apt install -yq bamtools && \
    apt clean all
	
USER ${NB_UID}

RUN mamba install --quiet -c bioconda -c anaconda --yes \
    biopython && \
    mamba clean  --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

WORKDIR "${HOME}"

USER root

# braker including RNAseq test file

RUN cd /opt && \
    git clone   --depth=1         https://github.com/Gaius-Augustus/BRAKER.git && \
    cd BRAKER && \
    cd example && \
    wget http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam

ENV PATH=${PATH}:/opt/BRAKER/scripts

# include ETP
RUN cd /opt && \
    git clone https://github.com/gatech-genemark/GeneMark-ETP.git && \
    mv GeneMark-ETP ETP && \
    chmod a+x /opt/ETP/bin/*py /opt/ETP/bin/*pl /opt/ETP/tools/*

ENV GENEMARK_PATH=/opt/ETP/bin
ENV PATH=${PATH}:/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes

USER ${NB_UID}

