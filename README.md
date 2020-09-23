# PipeRNAseq
A pipeline for RNAseq data analysis

## Software prerequisites
This pipeline is designed to run on Linux servers, and requires the following softwares:
```
bowtie2
STAR
bedtools
samtools
salmon
featureCount (from Subread)
fastqc (optional)
cufflinks (optional)
```
Besides the pipeline script PipeRNAseq.sh, dependencies are in ./bin folder

One UCSC tools (from http://hgdownload.soe.ucsc.edu/admin/exe/) is used: bedGraphToBigWig. Other scripts were generated from this project.

To save time, you can directly use STAR and featureCounts program in the ./bin folder (just add it to $PATH), without installing them again.

## Pipeline setup

Here is an example of mm10 genome setup. If you have set up PipeRiboseq.sh pipeline, then some folders don't need to set up again.

1, Download scripts from github to Linux server:

```
git clone https://github.com/sunyumail93/PipeRNAseq.git
mv PipeRNAseq PipelineHomeDir
```

2, Set up index files for genome mapping

2a, Download whole genome fasta sequence and chromosome sizes from UCSC goldenpath:

```
cd PipelineHomeDir/mm10/Sequence
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip *.gz
wget "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" -O mm10.ChromInfo.txt
```

2b, Extract RNA sequences
```
cd PipelineHomeDir/mm10/Annotation
gunzip *.gz
cd ../Sequence
bedtools getfasta -s -split -name -fi mm10.fa -bed ../Annotation/mm10.RefSeq.reduced.bed12 -fo mm10.RefSeq.reduced.bed12.fa.t
cat mm10.RefSeq.reduced.bed12.fa.t|sed 's/::.*//' > mm10.RefSeq.reduced.bed12.fa
rm -rf mm10.RefSeq.reduced.bed12.fa.t
```

2c, Set up index files:
```
mkdir PipelineHomeDir/mm10/Index
cd PipelineHomeDir/mm10/Index

#STAR index:
mkdir STARIndex
STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles ../Sequence/mm10.fa --sjdbGTFfile ../Annotation/mm10.RefSeq.reduced.bed12.geneid.gtf --sjdbOverhang 100

#salmon index (SalmonIndex directory will be created automatically):
#Genome FASTA index fai file will also be generated
salmon index -t ../Sequence/mm10.RefSeq.reduced.bed12.fa -i SalmonIndex --type quasi -k 31

#miRNA and rRNA bowtie2 index:
mkdir rRNAIndex
bowtie2-build ../Sequence/mm10.rRNA.fa ./rRNAIndex/rRNAIndex
```

4, Add executable permissions

```
chmod +x PipeRiboseq.sh
chmod +x ./bin/FastqAdapterTimmer
chmod +x ./bin/bedGraphToBigWig
```

## Pipeline components
```
PipelineHomeDir/
    ├── PipeRiboseq.sh
    ├── bin/
    └── mm10/
      └── Annotation/
        ├── mm10.RefSeq.reduced.bed12
        ├── mm10.RefSeq.reduced.mRNA.bed12
        ├── mm10.RefSeq.reduced.bed12.geneid.gtf
        └── mm10.uniqMatching.txt
      └── Index/
        ├── rRNAIndex/
        ├── SalmonIndex/
        └── STARIndex/
      └── Sequence/
        ├── mm10.fa
        ├── mm10.fai
        ├── mm10.ChromInfo.txt
        ├── mm10.rRNA.fa
        └── mm10.RefSeq.reduced.bed12.fa
    └── hg38/
       ...
```

Notes: 

1, For Annotation folder, download GTF file from UCSC table browser. `reduced`: Only one location was chosen when one gene duplicates at multiple genomic loci.

2, `uniqMatching.txt` file contains one-to-one matching from transcript to gene name.

3, For Index folder, indexes are not included in this github directory, but need to be created during set up.

4, For Sequence folder, `RefSeq.reduced.bed12.fa` was converted from `RefSeq.reduced.bed12` file using bedtools, and removed two rRNA genes (rRNA affects salmon quantification TPM). `genome.fa` and `ChromInfo.txt` files need to be downloaded from UCSC goldenpath. The fai index file is not required now sine it will be generated by samtools when needed.

## Usage

Type the pipeline name, then you will see the manual page:

```
PipeRNAseq.sh
```

Manual page:

![](images/Usages.png)

## Examples

A regular run using mostly default parameters:

```
#Single-end RNAseq
PipeRNAseq.sh -i Data.fastq.gz -g mm10
#Paired-end RNAseq
PipeRNAseq.sh -l Data.R1.fastq.gz -r Data.R2.fastq.gz -g mm10
```

More parameters used, run Cufflinks, not run fastqc, run unique mapping featureCouts (to pair with Riboseq data), generate bigWig tracks:

```
PipeRNAseq.sh -l Data.R1.fastq.gz -r Data.R2.fastq.gz -g -cufflinksrun -noqc -p 4 -bigWig -pairrpf
```

## Run a real data to test the pipeline

1, Download data

Use a public dataset: [GEO SRA: SRR10446759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4160756)

`fastq-dump` is part of [NCBI SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software):

```
fastq-dump --split-3 SRR10446759
gzip SRR10446759.fastq
```

2, Run PipeRNAseq.sh pipeline:

```
PipeRNAseq.sh -i SRR10446759.fastq.gz -g mm10
```

## Outputs

1, Gene expression
```
SRR10446759.fastq.mm10.featureCounts.gene.txt
SRR10446759.fastq.mm10.quant.sf
```

2, Tracks
```
SRR10446759.fastq.mm10.sorted.minus.bedGraph.bw
SRR10446759.fastq.mm10.sorted.minus.bedGraph.bw
```
