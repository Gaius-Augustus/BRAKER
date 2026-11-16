# Distributed under the terms of the Modified BSD License.
# Base: plain Ubuntu 22.04 (the jupyter/minimal-notebook image used before is
# also Ubuntu 22.04). Python and its packages come from a pixi environment
# (docker/pixi.toml, docker/pixi.lock) instead of conda.
ARG BASE_CONTAINER=ubuntu:22.04
FROM $BASE_CONTAINER as base

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    wget \
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

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

USER root

# basic tools that the jupyter base image used to provide
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    ca-certificates \
    git \
    wget \
    curl \
    bzip2 \
    xz-utils \
    zstd \
    unzip \
    procps \
    less && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=base /opt/ /opt/

ENV PATH=${PATH}:/opt/seqstats:/opt/cdbfasta:/opt/MakeHub

# AUGUSTUS does need several libraries that are now gone, re-install them:
RUN apt-get update --yes && \
    apt-get install -y libboost-iostreams-dev zlib1g-dev libboost-all-dev libboost-all-dev libbamtools-dev

ENV AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/

# augustus, install only in order to get the dependencies, will uninstall augustus later on
RUN apt update && \
    apt install -yq  augustus augustus-data augustus-doc \
    time && \
    apt clean all && \
    chmod -R a+rwX "${AUGUSTUS_CONFIG_PATH}"

# perl & dependencies
RUN apt update && \
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

# python and python packages: pixi environment from the lock file in
# /opt/pixi-env; /opt/conda links to it, so that scripts and workflows that
# put /opt/conda/bin on PATH (as BRAKER4 does) keep working.
ARG PIXI_VERSION=v0.81.0
ARG PIXI_SHA256=7aa3ec39aecceff9062fa2ed4d42cbaa0bdc25ddea727d048e061cf188d434f6
COPY docker/pixi.toml docker/pixi.lock /opt/pixi-env/
RUN cd /tmp && \
    wget -q https://github.com/prefix-dev/pixi/releases/download/${PIXI_VERSION}/pixi-x86_64-unknown-linux-musl.tar.gz && \
    echo "${PIXI_SHA256}  pixi-x86_64-unknown-linux-musl.tar.gz" | sha256sum -c - && \
    tar -xzf pixi-x86_64-unknown-linux-musl.tar.gz -C /usr/local/bin pixi && \
    rm pixi-x86_64-unknown-linux-musl.tar.gz && \
    cd /opt/pixi-env && \
    PIXI_CACHE_DIR=/tmp/pixi-cache pixi install --locked && \
    rm -r /tmp/pixi-cache && \
    ln -s /opt/pixi-env/.pixi/envs/default /opt/conda

ENV PATH=/opt/conda/bin:${PATH}

RUN apt-get remove -y augustus augustus-data augustus-doc

ENV AUGUSTUS_BIN_PATH=/opt/Augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/
ENV PATH=${PATH}:/opt/Augustus/scripts/:/opt/Augustus/bin/


# tsebra
RUN cd /opt && \
    git clone    https://github.com/Gaius-Augustus/TSEBRA

ENV PATH=${PATH}:/opt/TSEBRA/bin

# install genometools
RUN apt update && \
    apt install -yq genometools libgenometools0 libgenometools0-dev && \
    apt clean all

# makehub
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/MakeHub.git && \
    cd MakeHub && \
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
    chmod a+x bedToBigBed genePredCheck faToTwoBit gtfToGenePred hgGcPercent ixIxx  twoBitInfo wigToBigWig genePredToBed genePredToBigGenePred make_hub.py

ENV PATH=${PATH}:/opt/MakeHub

# perl dependencies of BRAKER and GeneMark-ETP+

RUN apt update && \
    apt install -yq libyaml-perl libhash-merge-perl libparallel-forkmanager-perl libscalar-util-numeric-perl libclass-data-inheritable-perl libexception-class-perl libtest-pod-perl libfile-which-perl libmce-perl libthread-queue-perl libmath-utils-perl libscalar-list-utils-perl && \
    apt clean all

# bamtools (ETP+)

RUN apt update && \
    apt install -yq bamtools && \
    apt clean all

# compleasm
RUN cd /opt && \
    wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.8/compleasm-0.2.8_x64-linux.tar.bz2 && \
    tar -xvjf compleasm-0.2.8_x64-linux.tar.bz2 && \
    rm compleasm-0.2.8_x64-linux.tar.bz2

# braker including RNAseq test file

RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/BRAKER.git && \
    cd BRAKER && \
    cd example && \
    wget http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam

ENV PATH=${PATH}:/opt/BRAKER/scripts

# include ETP (co-assembly fork, pinned to a commit)
ARG GENEMARK_ETP_COMMIT=8fcfe6f
RUN cd /opt && \
    git  clone https://github.com/Gaius-Augustus/GeneMark-ETP-co-assembly.git && \
    mv GeneMark-ETP-co-assembly GeneMark-ETP && \
    cd GeneMark-ETP && \
    git checkout ${GENEMARK_ETP_COMMIT} && \
    cd .. && \
    mv GeneMark-ETP ETP && \
    chmod a+x /opt/ETP/bin/*py /opt/ETP/bin/*pl /opt/ETP/tools/*

# gm-binaries-rust: unofficial Rust re-implementations of gmhmme3, gmhmmp and
# probuild (https://github.com/binaria-portata/gm-binaries-rust), CC BY-NC-SA
# 4.0, non-commercial use only, like GeneMark-ETP itself. They replace the
# GeneMark-ETP binaries under the original names; the originals are kept in
# /opt/gmhmm-orig (same layout as /opt/ETP/bin), the license and README of
# the re-implementation in /opt/gm-binaries-rust.
ARG GM_BINARIES_RUST_VERSION=v0.1.0
ARG GM_BINARIES_RUST_SHA256=1c9c7c56bef664182795b95388968c2a287b84772b4070ab17e4ac378c7dcebc
RUN cd /tmp && \
    d=gm-binaries-rust-${GM_BINARIES_RUST_VERSION}-x86_64-linux && \
    wget -q https://github.com/binaria-portata/gm-binaries-rust/releases/download/${GM_BINARIES_RUST_VERSION}/$d.tar.gz && \
    echo "${GM_BINARIES_RUST_SHA256}  $d.tar.gz" | sha256sum -c - && \
    tar -xzf $d.tar.gz && \
    cd /opt/ETP/bin && \
    mkdir -p /opt/gmhmm-orig /opt/gm-binaries-rust && \
    cp -p --parents gmes/gmhmme3 gmst/gmhmmp probuild gmes/probuild gmst/probuild gmes/ProtHint/dependencies/probuild /opt/gmhmm-orig/ && \
    install -m 755 /tmp/$d/gmhmm-rs gmes/gmhmme3 && \
    install -m 755 /tmp/$d/gmhmmp-rs gmst/gmhmmp && \
    for f in probuild gmes/probuild gmst/probuild gmes/ProtHint/dependencies/probuild; do install -m 755 /tmp/$d/probuild-rs $f; done && \
    cp /tmp/$d/LICENSE /tmp/$d/README.md /opt/gm-binaries-rust/ && \
    rm -r /tmp/$d /tmp/$d.tar.gz

# include EukSpecies-BRAKER2 because we use this container for braker-snake to prepare annotations
RUN cd /opt && \
    git clone https://github.com/gatech-genemark/EukSpecies-BRAKER2.git

ENV GENEMARK_PATH=/opt/ETP/bin
ENV PATH=${PATH}:/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes:/opt/compleasm_kit:/opt/EukSpecies-BRAKER2/bin
