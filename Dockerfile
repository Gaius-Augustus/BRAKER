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

# cdbfasta
RUN cd /opt && \ 
    git clone https://github.com/gpertea/cdbfasta.git && \
    cd cdbfasta && \
    make 
    
# get AUGUSTUS compilation dependencies
# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y samtools libhts-dev

# Install additional dependencies for homGeneMapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# compile augustus from source because of segmentation fault
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && \
    make && \
    cd scripts && \
    chmod a+x *.pl && \
    chmod a+x *.py


FROM $BASE_CONTAINER

USER root

COPY --from=base /opt/ /opt/

ENV PATH=${PATH}:/opt/seqstats:/opt/cdbfasta:/opt/MakeHub

# AUGUSTUS does need several libraries that are now gone, re-install them:
RUN apt-get update --yes && \
    apt-get install -y libboost-iostreams-dev zlib1g-dev libboost-all-dev libboost-all-dev libbamtools-dev

ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/

# augustus, install only in order to get the dependencies, will uninstall augustus later on
RUN apt update && \ 
    apt install -yq  augustus augustus-data augustus-doc \
    # for latex labels
    cm-super \
    dvipng \
    # for matplotlib anim
    ffmpeg \
    time && \
    apt clean all && \
    fix-permissions "${AUGUSTUS_CONFIG_PATH}"

# perl & dependencies
RUn apt update && \
    apt install -yq libyaml-perl \
                    libhash-merge-perl \
                    libparallel-forkmanager-perl \
                    libscalar-util-numeric-perl \
                    libclass-data-inheritable-perl \
                    libexception-class-perl \
                    libtest-pod-perl \
                    libfile-which-perl \
                    libmce-perl \
                    libthread-queue-perl \
                    libmath-utils-perl \
                    libscalar-list-utils-perl && \
    apt clean all

USER ${NB_UID}

# only python installations can be done as a normal user
RUN mamba install --quiet -c bioconda -c anaconda --yes \
    biopython && \
    mamba clean  --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

USER root

RUN apt-get remove -y augustus augustus-data augustus-doc

ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/
ENV PATH=${PATH}:/opt/Augustus/scripts/:/opt/Augustus/bin/


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

# perl dependencies of BRAKER and GeneMark-ETP+

RUn apt update && \
    apt install -yq libyaml-perl libhash-merge-perl libparallel-forkmanager-perl libscalar-util-numeric-perl libclass-data-inheritable-perl libexception-class-perl libtest-pod-perl libfile-which-perl libmce-perl libthread-queue-perl libmath-utils-perl libscalar-list-utils-perl && \
    apt clean all
    
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

RUN mamba install --quiet -c anaconda --yes \
    pandas && \
    mamba clean  --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

WORKDIR "${HOME}"

USER root

# compleasm
RUN cd /opt && \
    wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.5/compleasm-0.2.5_x64-linux.tar.bz2 && \
    tar -xvjf compleasm-0.2.5_x64-linux.tar.bz2 && \
    rm compleasm-0.2.5_x64-linux.tar.bz2

# braker including RNAseq test file

RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/BRAKER.git && \
    cd BRAKER && \
    cd example && \
    wget http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam

ENV PATH=${PATH}:/opt/BRAKER/scripts

# include ETP
RUN cd /opt && \
    git clone https://github.com/KatharinaHoff/GeneMark-ETP.git && \
#    cd GeneMark-ETP && \ # these lines are for the isoseq container
#    git checkout longread_experimental_dev && \
#    cd .. && \
    mv GeneMark-ETP ETP && \
    chmod a+x /opt/ETP/bin/*py /opt/ETP/bin/*pl /opt/ETP/tools/*

ENV GENEMARK_PATH=/opt/ETP/bin
ENV PATH=${PATH}:/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes:/opt/compleasm_kit

USER ${NB_UID}

