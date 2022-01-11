# Long-Read Integration Using BRAKER, GeneMarkS-T and TSEBRA

## Introduction

The following instructions are a prelimary protocol for producing a BRAKER<sup name="a1">[1](#ref1)</sup><sup name="a2">[2](#ref2)</sup><sup name="a3">[3](#ref3)</sup> gene set that integrates extrinsic evidence from <ins>long-read RNA-Seq</ins> (PacBio ccs), short-read RNA-seq (Illumina) and a large database of protein sequences (e.g. OrthoDB <sup name="a9">[9](#ref9)</sup> clade) into a single prediction. 

> :warning: This protocol was developed for presentation of intermediate results at **The Plant and Animal Genome Conference XXIX** (2022) in San Diego. Significant changes are expected prior publication in a peer reviewed scientific journal. For preliminary accuracy results, see https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/slides/slides_PAG2022.pdf (will be made available after January 9th 2022).

## Installation

#### BRAKER

You have to install the 'long_reads' branch of BRAKER via GitHub <https://github.com/Gaius-Augustus/BRAKER>. 

> :warning: The long read protocol does not work with e.g. an Anaconda installation of BRAKER because the long read branch is not part of the BRAKER release, yet.

Clone BRAKER:

```
git clone https://github.com/Gaius-Augustus/BRAKER
```

Switch to the 'long_reads' branch:

```
cd BRAKER/
git checkout long_reads
```

Follow the installation instructions of the README.md in the BRAKER repository.

Importantly, export ```BRAKER/scripts``` to your ```$PATH``` variable, e.g:

```
export PATH="/your/path/to/BRAKER/scripts:$PATH"
```
Make sure that ```braker.pl``` and ```gmst2globalCoords.py``` are executable (e.g. run ```which braker.pl``` and ```which gmst2globalCoords.py``` which should return the location of the respective scripts).

#### TSEBRA

You have to install the 'long_reads' branch of TSEBRA<sup name="a4">[4](#ref4)</sup> via GitHub <https://github.com/Gaius-Augustus/TSEBRA>.

> :warning: The long read protocol does not work with e.g. an Anaconda installation of TSEBRA because the long read branch is not part of the TSEBRA release, yet.

Clone TSEBRA:

```
git clone https://github.com/Gaius-Augustus/TSEBRA
```

Switch to 'long_reads' branch:

```
cd TSEBRA/
git checkout long_reads
```

Export ```TSEBRA/bin``` to your ```$PATH``` variable, e.g:

```
export PATH="/your/path/to/TSEBRA/bin:$PATH"
```

Make sure that ```tsebra.py``` is executable (e.g. by running ```which tsebra.py```, this should return the location of the tsebra.py executable).

#### GeneMarkS-T

Download GeneMarkS-T<sup name="a5">[5](#ref5)</sup> from <http://exon.gatech.edu/GeneMark/license_download.cgi> (the GeneMarkS-T) option. Unpack and install GeneMarkS-T as described in GeneMarkS-T’s `README` file.

Export ```gmst.pl``` to your ```$PATH``` variable, e.g:

```
export PATH="/your/path/to/gmst/gmst.pl:$PATH"
```

#### Minimap2
Download Minimap2<sup name="a6">[6](#ref6)</sup> from GitHub at <https://github.com/lh3/minimap2>. Follow the instructions in
the README file to install it and make sure that ```minimap2``` is executable.

#### Cupcake
You have to install Cupcake as it is described [here](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step#install).

Make sure that ```collapse_isoforms_by_sam.py``` is executable.

#### AUGUSTUS scripts

AUGUSTUS<sup name="a7">[7](#ref7)</sup><sup name="a8">[8](#ref8)</sup> should already be downloaded for the BRAKER installation. Make sure that the scripts at ```AUGUSTUS/scripts``` are executable (i.e. ```stringtie2fa.py```) and available in your ```$PATH``` variable, e.g:

```
export PATH="/your/path/to/AUGUSTUS/scripts:$PATH"
```

## Input Data

Below is a complete list of all necessary input files and their names used in this manual.

* ```genome.fa``` - genome sequence in FASTA format that has been softmasked for repeats
* ```RNAseq.bam``` - spliced alignments of short-read RNA-Seq in BAM format (RNA-Seq hints can be used instead of the BAM file, see the documentation of [BRAKER](https://github.com/Gaius-Augustus/BRAKER#data-description) for usage information)
* ```proteins.fa``` - a large database of protein sequences in FASTA format (e.g. a suitable OrthoDB partition)
* ```subreads1.fastq,subreads2.fastq,subreads3.fastq,...``` - list of **assembled** subread libraries from long-read RNA-Seq as FASTQ files (protocol has only been tested with PacBio ccs reads)

## Workflow

The workflow of this protocol is divided into four parts. In the first three steps, you have to generate three different gene sets, each using one source of extrinsic evidence. In the final step, these gene sets are combined using TSEBRA. The first three steps (BRAKER1, BRAKER2, GeneMarkS-T protocol) can be executed in any order.

Before you start, create a working directory for your project and set the number of threads that you want to use on your processor:

```
threads=8
wdir='/your/path/to/the/working/directory/'
mkdir $wdir
```

#### BRAKER1

Make sure that your ```$GENEMARK_PATH```(pointing to GeneMark-ES/ET/EP, not to GeneMarkS-T) and ```$AUGUSTUS_CONFIG_PATH``` paths are set correctly
or pass them manually to BRAKER by adding them in the command line below. See the BRAKER documentation for more information.

For convenience we will store all results in the same location, referred to as $wdir. You need to adapt this bash variable:

```
wdir=/a/directory/of/your/choice/
```

Run BRAKER using the spliced alignments of short-read RNA-seq to execute the BRAKER1 protocol:

```
braker.pl --genome=genome.fa --softmasking --cpu=$threads --bam=RNAseq.bam  --workingdir=$wdir/braker1/ 2> $wdir/braker1.log
```

The important results from this run are the BRAKER1 gene set ```augustus.hints.gtf``` and the hints from short-reads ```hintsfile.gff```. Make sure that they are available at ```$wdir/braker1/``` .

#### BRAKER2

Make sure that your ```$GENEMARK_PATH```, ```PROTHINT_PATH```, ```DIAMOND_PATH``` and ```$AUGUSTUS_CONFIG_PATH``` paths are set correctly
or pass them manually to BRAKER by adding them in the command line below. See the BRAKER documentation for more information.

Run BRAKER using the proteins to execute the BRAKER2 protocol:

```
braker.pl --genome=genome.fa --softmasking --cpu=$threads --epmode --prot_seq=proteins.fa  --workingdir=$wdir/braker2/ 2> $wdir/braker2.log
```

The important results from this run are the BRAKER2 gene set ```augustus.hints.gtf``` and the hints generated from the protein database ```hintsfile.gff```. Make sure that they are available at ```$wdir/braker2/``` .

#### GeneMarkS-T protocol

In this protocol, you map the transcripts from the assembled subread libraries from long-read RNA-Seq to the genome sequence and use GeneMarkS-T to identify protein-coding regions in the transcripts.

Create a working directory and change into it:

```
mkdir $wdir/long_read_protocol
cd $wdir/long_read_protocol
```

Concatenate the subread libraries files into a single file (note that you have to replace ```subreads1.fastq,subreads2.fastq,subreads3.fastq``` with your list of files in following command, it can be arbitrariely many files):

```
cat subreads1.fastq,subreads2.fastq,subreads3.fastq > all.subreads.fastq
```

Use Minimap2 to map the transcripts to the genome sequence, and then sort the result:

```
minimap2 -t $threads -ax splice:hq genome.fa all.subreads.fastq  > long_reads.sam
sort -k 3,3 -k 4,4n $wlong_reads.sam > long_reads.s.sam
```

Collapse redundant isoforms using a script from Cupcake:

```
collapse_isoforms_by_sam.py --input all.subreads.fastq --fq  -s long_reads.s.sam --dun-merge-5-shorter -o cupcake
```

You will find the collapsed transcripts in ```cupcake.collapsed.gff``` .

Run GeneMarkS-T to predict protein-coding regions in the transcripts:

```
stringtie2fa.py -g genome.fa -f cupcake.collapsed.gff -o cupcake.fa
gmst.pl --strand direct cupcake.fa.mrna --output gmst.out --format GFF
```

Use the GeneMarkS-T coordinates and the long-read transcripts to create a gene set in GTF format.

```
gmst2globalCoords.py -t cupcake.collapsed.gff -p gmst.out -o gmst.global.gtf -g genome.fa
```

The result of this step is located at ```$wdir/long_read_protocol/gmst.global.gtf``` .

#### TSEBRA

In the final step, the long-read version of TSEBRA is used to combine the three gene sets using all extrinsic evidence.

Create a working directory:

```
mkdir $wdir/tsebra
cd $wdir/tsebra
```

Run TSEBRA (note that you have to change the TSEBRA path to your TSEBRA installation directory for the ```-c``` option):

```
tsebra.py -g $wdir/braker1/augustus.hints.gtf,$wdir/braker2/augustus.hints.gtf -e $wdir/braker1/hintsfile.gff,$wdir/braker2/hintsfile.gff -l $wdir/long_read_protocol/gmst.global.gtf -c /your/path/to/TSEBRA/config/long_reads.cfg -o tsebra.gtf
```

The final gene set is located at ```$wdir/tsebra/tsebra.gtf``` .

## References

<b id="ref1">[1]</b> K. J. Hoff, S. Lange, A. Lomsadze, M. Borodovsky, and M. Stanke. 2015. “BRAKER1: Unsupervised Rna-Seq-Based Genome Annotation with GeneMark-Et and AUGUSTUS.” *Bioinformatics* 32 (5). Oxford University Press: 767--69.[↑](#a1)

<b id="ref2">[2]</b> T. Bruna, K. J. Hoff, A. Lomsadze, M. Stanke and M. Borodvsky. 2021. “BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database." *NAR Genomics and Bioinformatics* 3(1):lqaa108.[↑](#a2)

<b id="ref3">[3]</b> K. J. Hoff, A. Lomsadze, M. Borodovsky, and M. Stanke. (2019). Whole-Genome Annotation with BRAKER. Methods Mol Biol. 1962:65-95, doi: 10.1007/978-1-4939-9173-0_5. [↑](#a3)

<b id="ref4">[4]</b> L. Gabriel, K. J. Hoff, T. Bruna, M. Borodovsky and M. Stanke. 2021. "TSEBRA: transcript selector for BRAKER." *Bioinformatics* 22: 566.[↑](#a4)

<b id="ref5">[5]</b> S. Tang, A. Lomsadze and M. Borodovsky. 2015. "Identification of protein coding regions in RNA transcripts." *Nucleic acids research* 43 (12): e78-e78.[↑](#a5)

<b id="ref6">[6]</b>H. Li. 2018. "Minimap2: pairwise alignment for nucleotide sequences". *Bioinformatics* Volume 34 (18): 3094-3100.[↑](#a6)

<b id="ref7">[7]</b> M. Stanke, M. Diekhans, R. Baertsch, and D. Haussler. 2008. “Using Native and Syntenically Mapped cDNA Alignments to Improve de Novo Gene Finding.” *Bioinformatics* 24 (5). Oxford University Press: 637--44.[↑](#a7)

<b id="ref8">[8]</b> M. Stanke, O. Schöffmann, B. Morgenstern, and S. Waack. 2006. “Gene Prediction in Eukaryotes with a Generalized Hidden Markov Model That Uses Hints from External Sources.” *BMC Bioinformatics* 7 (1). BioMed Central: 62.[↑](#a8)

<b id="ref9">[9]</b> E. V. Kriventseva, D. Kuznetsov, F. Tegenfeldt, M. Manni, R. Dias, F. A. Simão and E. M. Zdobnov. 2019. OrthoDB v10: sampling the diversity of animal, plant, fungal, protist, bacterial and viral genomes for evolutionary and functional annotations of orthologs. *Nucleic acids research*, 47(D1), D807-D811.[↑](#a9)
