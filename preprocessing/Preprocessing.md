# Preprocessing of the annotation file
This protocol shows the procedures to prepare the annotation files, from UCSC gene annotation, using mm10 as an example.

All files and scripts used here are in: https://github.com/sunyumail93/PipeRNAseq/tree/master/preprocessing

## Download annotation files
Download annotation from UCSC Genome Browser, [Table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1269791857_CW9VuvTqCYCCUG0WVAGccGPx7DbS).

Select:

group: Genes and Gene Predictions

track: NCBI RefSeq

table: RefSeq Curated (ncbiRefSeqCurated)

Then download 'all fields from selected table' as `mm10.RefSeq.UCSC.all`

Also download 'GTF - gene transfer format (limited)' as `mm10.RefSeq.UCSC.gtf`

This GTF file can be renamed as `mm10.RefSeq.simplified.geneid.gtf` and used in PipeRNAseq.sh. However, this GTF file contains duplicated transcripts, for example, transcript `NR_165495` has 10 copies on chr7, and transcript `NR_162797` has 3 copies on differnt chromosomes:

```
$ grep NR_165495 mm10.RefSeq.UCSC.gtf
chr7	mm10_refGene	exon	59466292	59466370	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495"; 
chr7	mm10_refGene	exon	59462548	59462626	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup1"; 
chr7	mm10_refGene	exon	59464420	59464498	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup2"; 
chr7	mm10_refGene	exon	59468163	59468241	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup3"; 
chr7	mm10_refGene	exon	59393561	59393639	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup4"; 
chr7	mm10_refGene	exon	59412802	59412880	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup5"; 
chr7	mm10_refGene	exon	59415590	59415668	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup6"; 
chr7	mm10_refGene	exon	59417457	59417535	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup7"; 
chr7	mm10_refGene	exon	59421187	59421265	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup8"; 
chr7	mm10_refGene	exon	59460677	59460755	0.000000	-	.	gene_id "NR_165495"; transcript_id "NR_165495_dup9"; 

$ grep NR_162797 mm10.RefSeq.UCSC.gtf
chr2	mm10_refGene	exon	180894045	180894101	0.000000	-	.	gene_id "NR_162797"; transcript_id "NR_162797"; 
chr3	mm10_refGene	exon	17795685	17795741	0.000000	-	.	gene_id "NR_162797"; transcript_id "NR_162797"; 
chr14	mm10_refGene	exon	64590669	64590727	0.000000	-	.	gene_id "NR_162797"; transcript_id "NR_162797"; 
```

We can count how many duplications there:

```
awk 'NR>1' mm10.RefSeq.UCSC.all|wc -l                                 47382       (All transcripts, including duplication)
awk '{print $10}' mm10.RefSeq.UCSC.gtf|sort|uniq|wc -l                45895       Transcripts (transcript duplications removed)
awk 'NR>1' mm10.RefSeq.UCSC.all|awk '{print $13}'|sort|uniq|wc -l     26214       Genes
```

This file will cause error when converting GTF to bed12 format using `gtfToGenePred` and `genePredToBed` from [UCSC tools](http://hgdownload.soe.ucsc.edu/admin/exe/).

So I used the following script to remove duplications and only keep one copy of each transcript:

```
$ chmod +x RefSeqGTF_Simplifier.sh
$ RefSeqGTF_Simplifier.sh
[Usage]: RefSeqGTF_Simplifier.sh [Genome.RefSeq.gtf]

# This command may run ~15 hours, and it will generate mm10.RefSeq.UCSC.simplified.gtf
RefSeqGTF_Simplifier.sh mm10.RefSeq.UCSC.gtf
```

Also prepare the uniqMatching file to match transcript ID and gene name:

```
awk 'NR>1' mm10.RefSeq.UCSC.all |awk 'BEGIN{OFS="\t"}{print $2,$13}'|sort|uniq > mm10.uniqMatching_WithrRNA.txt
wc -l mm10.uniqMatching_WithrRNA.txt                                           45895
awk '{print $2}' mm10.uniqMatching_WithrRNA.txt|sort|uniq > mm10.uniqMatching.gene_WithrRNA.txt
```

Add the Gene name into the GTF file:
```
chmod +x RefSeqGTF_NameAdder.sh
$ RefSeqGTF_NameAdder.sh
Usage: RefSeqGTF_NameAdder.sh [Ori.gtf] [Unique Matching|oriName newName] [New.gtf]

RefSeqGTF_NameAdder.sh mm10.RefSeq.UCSC.simplified.gtf mm10.uniqMatching_WithrRNA.txt mm10.RefSeq.reduced.bed12.geneid_WithrRNA.gtf
```

Finally, remove the rRNA genes, since usually due to poly-A selection or rRNA depletion method, rRNA genes cannot be quantified accurately.

Different species may have different names of rRNA genes. In mm10, they are Rn4.5s and Rn45s:

```
# There are two rRNA genes in mm10 annotation:
NR_002841	Rn4.5s
NR_046233	Rn45s

grep -v "NR_002841\|NR_046233" mm10.RefSeq.reduced.bed12.geneid_WithrRNA.gtf > mm10.RefSeq.reduced.bed12.geneid.gtf
grep -v "NR_002841\|NR_046233" mm10.uniqMatching_WithrRNA.txt > mm10.uniqMatching.txt

#Use gtfToGenePred from http://hgdownload.soe.ucsc.edu/admin/exe/
gtfToGenePred mm10.RefSeq.reduced.bed12.geneid.gtf mm10.RefSeq.reduced.bed12.geneid.gp
genePredToBed mm10.RefSeq.reduced.bed12.geneid.gp mm10.RefSeq.reduced.bed12

#Extract mRNA NM sequences:
awk '$4~/^NM/' mm10.RefSeq.reduced.bed12 > mm10.RefSeq.reduced.mRNA.bed12
```

If you need to extract CDS, you can use the BED12Extractor.sh here: https://github.com/sunyumail93/Bed12Processing

```
BED12Extractor.sh -a cds -i mm10.RefSeq.reduced.mRNA.bed12 -o mm10.RefSeq.reduced.mRNA.cds.bed12
```
