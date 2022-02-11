# pRNASeqTools
Integrated High-throughput Sequencing Data Analysis for Plant

Author: Dr. Chenjiang You

Current Version: 0.8

Latest updata: 02/09/2022

- - - -
### Introduction

The `pRNASeqTools` is a `Perl` and `R` based pipeline designed for automatic general analysis for _Illumina_ sequencing data in supported plants.

Currently it is able to process small RNA-seq, mRNA-seq, degradome-seq, CLIP-seq, ChIP-seq, and WGBS-seq and generally analyze them. See below for more information of specific tasks.

The phasiRNA identification module was contributed by [Dr. Xuan Ma](mailto:xuanma@genetics.ac.cn).

If you have any questions or comments, please submit an issue in the GitHub or directly email to [Chenjiang You](mailto:cjyou@fudan.edu.cn).

- - - -
### Before the first run

To successfully run this pipeline on your own computer or server, several pieces of dependent software are needed. See `INSTALL.md` for detailed information.

For genome reference files, please contact [Chenjiang You](mailto:cjyou@fudan.edu.cn) for pre-built genomes or the instructions for new genomes.

- - - -
### Input files preparation

The only input files needed are _Illumina_ output fastq files, either in the `FASTQ` format or corresponding compressed file formats `.gz` and `.bz2`. SRR accessions are also accepted.
**Note: the `FASTA` format is not supported. You may convert the fasta format to fastq format by adding artificial sequence names and qualities.**

- - - -
### Running pRNASeqTools

#### Usage and Options

See help information of `pRNASeqTools` simply by execute `pRNASeqTools`.

### Quick examples

General analysis for small RNA-seq from samples `control` and `treatment` with 3 biological replicates

```bash
pRNASeqTools srna --adaptor AGATCGGAAGAGC --control control=control_1.fastq.gz+control_2.fastq.gz+control_3.fastq.gz --treatment treatment=treatment_1.fastq.bz2+treatment_2.fastq.bz2+treatment_3.fastq.bz2
```

Only mapping the small RNA reads to the genome and creating read count files

```bash
pRNASeqTools srna --adaptor AGATCGGAAGAGC --mapping-only --control control=control_1.fastq+control_2.fastq+SRRXXXXXXX
```

Perform statistic analyses in the folder containing pre-processed data

```bash
pRNASeqTools srna --nomapping --control control=3 --treatment treatment=3
```

General analysis for mRNA-seq from samples `control` and `treatment` with 3 biological replicates

```bash
pRNASeqTools mrna --control control=control_1.fastq.gz+control_2.fastq.gz+control_3.fastq.gz --treatment treatment=treatment_1.fastq.bz2+treatment_2.fastq.bz2+treatment_3.fastq.bz2
```

General analysis for paired mRNA-seq from samples `control` and `treatment` with 3 biological replicates

```bash
pRNASeqTools mrna --control control=control_1_R1.fastq.gz,control_1_R2.fastq.gz+control_2_R1.fastq.gz,control_2_R2.fastq.gz+control_3_R1.fastq.gz,control_3_R2.fastq.gz --treatment treatment=treatment_1.fastq.bz2+treatment_2.fastq.bz2+treatment_3.fastq.bz2
```

Trunction and tailing analysis of plant miRNAs

```bash
pRNASeqTools tt --adaptor AGATCGGAAGAGC --control control=control_1.fastq.gz+control_2.fastq.gz+control_3.fastq.gz --treatment treatment=treatment_1.fastq.bz2+treatment_2.fastq.bz2+treatment_3.fastq.bz2
```

Degradome data analysis

```bash
pRNASeqTools deg --adaptor AGATCGGAAGAGC --control control=control_1.fastq.gz+control_2.fastq.gz+control_3.fastq.gz --treatment treatment=treatment_1.fastq.bz2+treatment_2.fastq.bz2+treatment_3.fastq.bz2
```

Two-factor DE analysis

```bash
pRNASeqTools tf --control control=time1,3,time2,3 --treatment treatment=time1,3,time2,3
```

- - - -
### Output

All output files are stored in the output directory.
Mapping statistics are stored in the log file `log_xxxxxxxxx.txt`.

#### sRNA analysis

Several groups of files are generated in the output directory:

1. `count` files and `nf` files can be used for later `--nomapping` runs, which will not invoke the mapping procedures.

   **The second to tenth columns of `count` files are numbers of assigned small RNAs with length 18 - 26nt.**

2. `pdf` files showing the reproductivity of biological replicates and the relationship of samples.

3. `csv` files containing the results of statistic analyses, of which `hyper` and `hypo` files indicate the significant ones filtered out based on input parameters.

4. `bedgraph` files for visualization in IGV.
   **Note: Keywords are embedded in the file names, indicating the targets and methods.**

#### phased siRNA analysis

1.

#### Trunction and tailing analysis

1. miRNA reads are categorized in the `out` files.
__The second column shows the number of tailed nucleotides and the third column shows the number of truncated nucleotides.__
2. `pdf` files are bubble plots for each miRNA.

#### Degradome data analysis

1. `bam` files contain the mapped reads.
2. `txt` files report the identified peaks on each transcripts.

#### mRNA analysis

1. Mapped reads in `bam` files and read counts for each gene in `txt` files are reported.
2. Up-regulated and down-regulated DEG results are reported in `total.hyper.csv` and `total.hypo.csv` files.

#### Ribo-seq analysis

#### Two-factor of DEG analysis

1. This mode can only run in `srna` and `mrna` output folders.
2. Up-regulated and down-regulated DEG results are reported in `total.hyper.csv` and `total.hypo.csv` files.

#### CLIP-seq analysis

#### ChIP-seq analysis

#### WGBS-seq analysis

#### RiboMeth-seq analysis
