#!/bin/bash
#PipeRNAseq.sh
#This is a pipeline for RNAseq FastQ analysis, version 6
#Inputs: RNAseq raw reads in FASTQ format, also the genome to be used in analysis
#Outputs: 
#Usage1: piRNA_Targets.sh -i [Data.fastq] -g mm10
#Usage2: piRNA_Targets.sh -l [Data.R1.fastq] -r [Data.R2.fastq] -g hg38
#Example: PipeRNAseq.sh -i Embryo2cell.fastq.gz -g mm10
#Example: PipeRNAseq.sh -l B6.Testis.R1.fastq -r B6.Testis.R2.fastq -g mm10
#Example: PipeRNAseq.sh -l B6.Testis.R1.fastq.gz -r B6.Testis.R2.fastq.gz -g mm10

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                       PipeRNAseq: pipeline for RNAseq analysis.                          *"
  echo "*                             Version 6, 2019-11-05, Y.S                                   *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (single end data):                                                       *"
  echo "*                  -i [Data.fastq]                                                         *"
  echo "*                  -g [mm10/hg38/rn6/GRCm38/SacCer3/susScr11]                              *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -l [Data.R1.fastq]                                                      *"
  echo "*                  -r [Data.R2.fastq]                                                      *"
  echo "*                  -g [mm10/hg38/rn6/GRCm38/SacCer3/susScr11]                              *"
  echo "*        Optional: -p [Number of CPUs, default=1]                                          *"
  echo "*                  -t [Cufflinks library type, default=fr-firststrand]                     *"
  echo "*                  -s [Salmon strand type, default=A, automatically detect]                *"
  echo "*                  -pre Sequences.fa [Run pre-mapping to Sequences.fa]                     *"
  echo "*                  -cufflinksrun [Run cufflinks]                                           *"
  echo "*                  -noqc [Suppress fastqc]                                                 *"
  echo "*                  -pairrpf [Run unique CDS featureCounts]                                 *"
  echo "*                  -bigWig [Generate bigWig by default, genomic mapping reads as norm]     *"
  echo "* Inputs: All data need to have .fastq (or.fastq.gz) as suffix for correct recognition     *"
  echo "* Run: Default to run fastqc, rRNA & genomic mapping, featureCounts and salmon             *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "********************************************************************************************"
  exit 1
fi

echo "*****************************************************************************"
echo "*                PipeRNAseq: pipeline for RNAseq analysis.                  *"
echo "*                      Version 6, 2019-11-05,  Y.S                          *"
echo "*****************************************************************************"
echo "0. Loading softwares:"
#Get current time
St=`date`
calculating(){ awk "BEGIN { print "$*" }"; }

#Load softwares: This part may need to be changed when using different cluster
#Main tools:bowtie2, STAR, cufflinks, salmon, featureCounts (from Subread package)
#Other tools: samtools, bedtools, fastqc, Python2

#Get pipeline directory
HomeDir=$(dirname `readlink -f $0`)
#For mac local check only:
#HomeDir="/Users/yusun/Downloads/PipelineHomeDir"
echo "   Home Directory:"
echo "   "$HomeDir

echo "1. Resolving inputs:"

#Get parameters
Data="unassigned"
DataLeft="unassigned"
DataRight="unassigned"
datastats=0
genome="unassigned"
genomestats=0
CPU=1
LibraryType="fr-firststrand"
StrandType="unassigned"
runfastqc=1              #You can change this default to allow/suppress fastqc
runcufflinks=0
premap=0
premapData="unassigned"
pairrpf=0
bigWig=1                 #Generate bigWig as default

for arg in "$@"
do
 if [[ $arg == "-i" ]]
  then
    Data=${array[$counter+1]}
    echo '   Single end data: '$Data
 elif [[ $arg == "-l" ]]    
  then
    DataLeft=${array[$counter+1]}
    echo '   Paired end data, left end: '$DataLeft
 elif [[ $arg == "-r" ]]
  then
    DataRight=${array[$counter+1]}
    echo '   Paired end data, right end: '$DataRight
 elif [[ $arg == "-g" ]]
  then
    genome=${array[$counter+1]}
    echo '   Genome: '$genome
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-t" ]]
  then
    LibraryType=${array[$counter+1]}
    echo '   Library type: '$LibraryType
 elif [[ $arg == "-s" ]]
  then
    StrandType=${array[$counter+1]}
    echo '   Salmon strand type: '$StrandType
 elif [[ $arg == "-pre" ]]
  then
    premapData=${array[$counter+1]}
    premap=1
    echo '   Pre mapping file: '$premapData
 elif [[ $arg == "-cufflinksrun" ]]
  then
    runcufflinks=1
    echo '   Run cufflinks'
 elif [[ $arg == "-noqc" ]]
  then
    runfastqc=0
    echo '   Suppress fastqc quality control'
 elif [[ $arg == "-pairrpf" ]]
  then
    pairrpf=1
    echo '   Run unique CDS featureCounts'
 elif [[ $arg == "-bigWig" ]]
  then
    bigWig=1
    echo '   Generate bigWig tracks, which requires bedGraphToBigWig'
 fi
  let counter=$counter+1
done

#Get current directory and create folders or files
[ $runfastqc == "1" ] && FastqcDir=fastqc && mkdir -p $FastqcDir
[ $premap == "1" ] && PremapDir=pre_mapping && mkdir -p $PremapDir
GenomeMappingDir=genome_mapping && mkdir -p $GenomeMappingDir
[ $runcufflinks == "1" ] && CufflinksDir=cufflinks_results && mkdir -p $CufflinksDir
SalmonOutputDir=salmon_results && mkdir -p $SalmonOutputDir
FeatureDir=feature_counts && mkdir -p $FeatureDir
[ -d $HomeDir/$genome/Index/SalmonIndexWithTE ] && SalmonOutputTEDir=salmonWithTE_results && mkdir -p $SalmonOutputTEDir
[ $bigWig == "1" ] && TracksDir=tracks && mkdir -p $TracksDir

#Check data and determine single or paired mode
if [ $Data == "unassigned" -a $DataLeft == "unassigned" -a $DataRight ==  "unassigned" ];then
  echo "     >>> [Error]: Please input data!"
else
  if [ $Data != "unassigned" ];then
    if [ $DataLeft == "unassigned" -a $DataRight ==  "unassigned" ];then
      echo "     >>> Run single end mode"
      mode="single"
      datastats=1
    else
      echo "     >>> [Error]: Mixed single and paired end inputs!"
    fi
  elif [ $DataLeft != "unassigned" -a $DataRight !=  "unassigned" ];then
    echo "     >>> Run paired end mode"
    mode="paired"
    datastats=1
  else
    echo "     >>> [Error]: Data not paired!"
  fi
fi

#Getting Output Suffix and default strand types
if [ $mode == "single" ];then
  OutputSuffix=$(echo $Data|sed 's/.fastq.*//g')
  TABLE=${OutputSuffix}.summary
  if [ $StrandType == "unassigned" ];then
    StrandType="A"
  fi
else
  OutputSuffix=$(echo $DataLeft|sed 's/.R[0-9].*fastq.gz$//g'|sed 's/.R[0-9].*fastq$//g')
  TABLE=${OutputSuffix}.summary
  if [ $StrandType == "unassigned" ];then
    StrandType="A"
  fi
fi

#Check genome
if [ $genome == "unassigned" ];then
  echo "     >>> [Error]: Please assign genome file!"
else
  if [ -d $HomeDir/$genome ];  then
    echo "     >>> This genome is supported."
    genomestats=1
  else
    echo "     >>> [Error]: Genome unsupported!"
  fi
fi

#Export checking status
if [ $datastats == "1" -a $genomestats ==  "1" ];then
    echo "   Output Suffix: "$OutputSuffix
  echo "   Status check: pass"
else
  echo "   Status check: failed, stop the pipeline"
  exit 1
fi

echo "*****************************************************************************"
echo "2. Checking dependencies:"
echo "Annotations:"
#Genome files and index check
#This check includes folders: /Annotation, /Sequence, /Index
#                    files: ${genome}.RefSeq.gtf, ${genome}.RefSeq.bed12
#                    files: $genome.fa, $genome.RefSeq.fa
dependenciescount=0

#Annotation and FASTA files
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12 ];then
  echo "   BED12 Annotation: "${genome}.RefSeq.reduced.bed12
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] BED12 Annotation lost "
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf ];then
  echo "   GTF Annotation: "${genome}.RefSeq.reduced.bed12.geneid.gtf
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] GTF Annotation lost "
fi
if [ -s $HomeDir/$genome/Annotation/${genome}.uniqMatching.txt ];then
  echo "   RefSeq to refFlat: "${genome}.uniqMatching.txt
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] RefSeq to refFlat Matching lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.fa ];then
  echo "   Genome: "${genome}.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Genome lost "
fi
if [ -s $HomeDir/$genome/Sequence/${genome}.RefSeq.reduced.bed12.fa ];then
  echo "   Transcriptome: "${genome}.RefSeq.reduced.bed12.fa
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] Transcriptome lost "
fi

###Index folders
echo "Index files:"
if [ -d $HomeDir/$genome/Index/STARIndex ];then
  echo "   STARIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] STARIndex lost"
fi
if [ -d $HomeDir/$genome/Index/SalmonIndex ];then
  echo "   SalmonIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] SalmonIndex lost "
fi
if [ -d $HomeDir/$genome/Index/rRNAIndex ];then
  echo "   rRNAIndex exists"
  let dependenciescount=$dependenciescount+1
else
  echo "   [Error] rRNAIndex lost "
fi

#echo $dependenciescount
if [ $dependenciescount == "8" ];then
  echo "   Dependencies check: pass"
else
  echo "   [Error] Some dependencies lost "
fi

echo "*****************************************************************************"
echo "3. Run fastqc quality control:"
if [ $runfastqc == "1" -a $mode == "single" ];then
echo "   Running fastqc single-end mode"
fastqc \
    -f fastq \
    -o $FastqcDir \
    $Data \
    2> $FastqcDir/${OutputSuffix}.fastqc.log && \
    echo "   Done fastqc"
elif [ $runfastqc == "1" -a $mode == "paired" ];then
echo "   Running fastqc paired-end mode"
fastqc \
    -f fastq \
    -o $FastqcDir \
    $DataLeft $DataRight \
    2> $FastqcDir/${OutputSuffix}.fastqc.log && \
    echo "   Done fastqc"
elif [ $runfastqc == "0" ];then
echo "   Skipping fastqc..."
fi

echo "*****************************************************************************"
echo "4. rRNA mapping for quality check:"
if [ $mode == "single" ];then
echo "   Running bowtie2 single-end mode"
bowtie2 \
    -x $HomeDir/$genome/Index/rRNAIndex/rRNAIndex \
    -U $Data \
    -q \
    --very-fast \
    -k 1 \
    -p $CPU \
    -S /dev/null \
    2> $GenomeMappingDir/${OutputSuffix}.rRNA.log && \
    echo "   Done rRNA mapping"
	rRNAReads=`head -4 $GenomeMappingDir/${OutputSuffix}.rRNA.log | tail -1 | awk '{print $1}'`
else
echo "   Running bowtie2 paired-end mode"
bowtie2 \
    -x $HomeDir/$genome/Index/rRNAIndex/rRNAIndex \
    -1 $DataLeft \
    -2 $DataRight \
    -q \
    --very-fast \
    -k 1 \
    --no-mixed \
    --no-discordant \
    -p $CPU \
    -S /dev/null \
    2> $GenomeMappingDir/${OutputSuffix}.rRNA.log && \
    echo "   Done rRNA mapping"
	rRNAReads=`head -4 $GenomeMappingDir/${OutputSuffix}.rRNA.log | tail -1 | awk '{print $1}'`
fi

if [ $premap == "1" ];then
echo "4.a Running bowtie2 single-end mode for extra sequences"
echo "   Building index for input sequence using bowtie2"
bowtie2-build $premapData $PremapDir/ExtraSeq 2> $PremapDir/${OutputSuffix}.index.log 1> $PremapDir/temp && rm -rf $PremapDir/temp
if [ $mode == "single" ];then
bowtie2 \
    -x $PremapDir/ExtraSeq \
    -U $Data \
    -q \
    -a \
    --very-fast \
    -p $CPU \
    -S $PremapDir/${OutputSuffix}.premap.sam \
    2> $PremapDir/${OutputSuffix}.premap.log && \
    echo "   Done extra sequences mapping"
	premapReads=`head -4 $PremapDir/${OutputSuffix}.premap.log | tail -1 | awk '{print $1}'`
else
bowtie2 \
    -x $PremapDir/ExtraSeq \
    -1 $DataLeft \
    -2 $DataRight \
    -q \
    -a \
    --very-fast \
    --no-mixed \
    --no-discordant \
    -p $CPU \
    -S $PremapDir/${OutputSuffix}.premap.sam \
    2> $PremapDir/${OutputSuffix}.premap.log && \
    echo "   Done extra sequences mapping"
    premapReads=`head -4 $PremapDir/${OutputSuffix}.premap.log | tail -1 | awk '{print $1}'`
fi
echo "   Summarizing extra sequences mapping results"
samtools view -F 4 $PremapDir/${OutputSuffix}.premap.sam > $PremapDir/${OutputSuffix}.premap.sam.hits
samtools view -bS $PremapDir/${OutputSuffix}.premap.sam > $PremapDir/${OutputSuffix}.premap.bam && rm -rf $PremapDir/${OutputSuffix}.premap.sam
rm -rf $PremapDir/${OutputSuffix}.premap.summary.txt && echo -e "Gene\tReads" > $PremapDir/${OutputSuffix}.premap.summary.txt
grep '>' $premapData | sed 's/>//' > $PremapDir/FastaList
awk '{print $3}' $PremapDir/${OutputSuffix}.premap.sam.hits |sort|uniq -c|awk '{OFS="\t";print $2,$1}' > $PremapDir/${OutputSuffix}.premap.summary.temp
awk '{print $1}' $PremapDir/${OutputSuffix}.premap.summary.temp > $PremapDir/${OutputSuffix}.premap.summary.tempList
cat $PremapDir/FastaList $PremapDir/${OutputSuffix}.premap.summary.tempList | sort|uniq -c|awk '$1==1'|awk '{OFS="\t";print $2,0}' > $PremapDir/${OutputSuffix}.premap.summary.tempfill
cat $PremapDir/${OutputSuffix}.premap.summary.temp $PremapDir/${OutputSuffix}.premap.summary.tempfill > $PremapDir/${OutputSuffix}.premap.summary.tempall
for name in `cat $PremapDir/FastaList`;do awk -v t=$name '{OFS="\t";if ($1==t) print $0}' $PremapDir/${OutputSuffix}.premap.summary.tempall;done >> $PremapDir/${OutputSuffix}.premap.summary.txt
rm -rf $PremapDir/${OutputSuffix}.premap.summary.temp* $PremapDir/FastaList 
fi

echo "*****************************************************************************"
echo "5. Genomic mapping using STAR:"
if [ $mode == "single" ];then
echo "   Running STAR single-end mode"
  if [[ $Data == *.gz ]];then
  echo "   Running zipped data"
STAR \
    --runMode alignReads \
    --genomeDir $HomeDir/$genome/Index/STARIndex \
    --readFilesIn $Data \
    --runThreadN $CPU \
    --readFilesCommand zcat \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMunmapped Within KeepPairs \
    --outSAMattributes All \
    --outFileNamePrefix $GenomeMappingDir/${OutputSuffix}.${genome}. \
    --outSAMtype BAM Unsorted 2>&1 1> $GenomeMappingDir/${OutputSuffix}.STAR.log && \
    echo "   Done STAR mapping"
  else
STAR \
    --runMode alignReads \
    --genomeDir $HomeDir/$genome/Index/STARIndex \
    --readFilesIn $Data \
    --runThreadN $CPU \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMunmapped Within KeepPairs \
    --outSAMattributes All \
    --outFileNamePrefix $GenomeMappingDir/${OutputSuffix}.${genome}. \
    --outSAMtype BAM Unsorted 2>&1 1> $GenomeMappingDir/${OutputSuffix}.STAR.log && \
    echo "   Done STAR mapping"
  fi
else
echo "   Running STAR paired-end mode"
  if [[ $DataLeft == *.gz ]];then
  echo "   Running zipped data"
STAR \
    --runMode alignReads \
    --genomeDir $HomeDir/$genome/Index/STARIndex \
    --readFilesIn $DataLeft $DataRight \
    --runThreadN $CPU \
    --readFilesCommand zcat \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMunmapped Within KeepPairs \
    --outSAMattributes All \
    --outFileNamePrefix $GenomeMappingDir/$OutputSuffix.${genome}. \
    --outSAMtype BAM Unsorted 2>&1 1> $GenomeMappingDir/${OutputSuffix}.STAR.log && \
    echo "   Done STAR mapping"
  else
STAR \
    --runMode alignReads \
    --genomeDir $HomeDir/$genome/Index/STARIndex \
    --readFilesIn $DataLeft $DataRight \
    --runThreadN $CPU \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outSAMunmapped Within KeepPairs \
    --outSAMattributes All \
    --outFileNamePrefix $GenomeMappingDir/$OutputSuffix.${genome}. \
    --outSAMtype BAM Unsorted 2>&1 1> $GenomeMappingDir/${OutputSuffix}.STAR.log && \
    echo "   Done STAR mapping"
  fi
fi

echo "*****************************************************************************"
echo "6. Post-mapping processing:"
echo "   Sorting the bam file..."
samtools sort -T $GenomeMappingDir/aln.sorted $GenomeMappingDir/${OutputSuffix}.${genome}.Aligned.out.bam -o $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam
echo "   Indexing..."
samtools index $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam
echo "   Getting unique mapping bam file..."
samtools view $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam
samtools view -H $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header
samtools view -q 10 $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10
cat $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10 > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.sam
samtools view -b $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.sam > $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.header $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam.uMAPQ10
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.sam $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.sam
#This un-sorted bam is not quite useful now
rm -rf $GenomeMappingDir/${OutputSuffix}.${genome}.Aligned.out.bam
samtools index $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam

echo "   Getting statistics:"
rm -rf $TABLE
InputReads=`grep 'Number of input reads' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
UniquReads=`grep 'Uniquely mapped reads number' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
MultiReads=`grep 'Number of reads mapped to multiple loci' $GenomeMappingDir/$OutputSuffix.${genome}.Log.final.out | awk '{print $NF}'`
AllMapReads=$((UniquReads+MultiReads))
UnMapReads=$((InputReads-UniquReads-MultiReads))
MappingRate=`awk -v mapping=$AllMapReads -v total=$InputReads 'BEGIN{print (mapping/total*100)}'`
ContaminationRate=`awk -v mapping=$rRNAReads -v total=$InputReads 'BEGIN{print (mapping/total*100)}'`
echo -e "   total input reads:\t${InputReads}"
echo -e "   rRNA_reads:\t${rRNAReads}"
if [ $premap == "1" ];then
  echo -e "   premap_reads:\t${premapReads}"
fi
echo -e "   genomic_mapped_reads:\t${AllMapReads}"
echo -e "   genomic_unique_mapped_reads:\t${UniquReads}"
echo -e "   genomic_multiple_mapped_reads:\t${MultiReads}"
echo -e "   genomic_unmappable_reads:\t${UnMapReads}"
echo -e "   rRNA mapping rate(%):\t${ContaminationRate}"
echo -e "   Mapping rate(%):\t${MappingRate}"
echo -e "total input reads:\t${InputReads}" >> $TABLE
echo -e "rRNA_reads:\t${rRNAReads}" >> $TABLE
if [ $premap == "1" ];then
  echo -e "premap_reads:\t${premapReads}" >> $TABLE
fi
echo -e "genomic_mapped_reads:\t${AllMapReads}" >> $TABLE
echo -e "genomic_unique_mapped_reads:\t${UniquReads}" >> $TABLE
echo -e "genomic_multiple_mapped_reads:\t${MultiReads}" >> $TABLE
echo -e "genomic_unmappable_reads:\t${UnMapReads}" >> $TABLE
echo -e "rRNA mapping rate(%):\t${ContaminationRate}" >> $TABLE
echo -e "Mapping rate(%):\t${MappingRate}" >> $TABLE

echo "*****************************************************************************"
echo "7. Transcriptome assembly and abundance calculation using cufflinks:"
if [ $runcufflinks == "1" -a $mode == "single" ];then
echo "   Running cufflinks single-end mode"
echo "   Using library type: "$LibraryType
cufflinks \
    -o $CufflinksDir \
    -p $CPU \
    -G $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -b $HomeDir/$genome/Sequence/${genome}.fa \
    -u \
    --library-type $LibraryType \
    --compatible-hits-norm \
    --no-update-check \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
                2> $CufflinksDir/${OutputSuffix}.${genome}.cufflinks.log
elif [ $runcufflinks == "1" -a $mode == "paired" ];then
echo "   Running cufflinks paired-end mode"
echo "   Using library type: "$LibraryType
cufflinks \
    -o $CufflinksDir \
    -p $CPU \
    -G $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -b $HomeDir/$genome/Sequence/${genome}.fa \
    -u \
    --library-type $LibraryType \
    --compatible-hits-norm \
    --no-update-check \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
                2> $CufflinksDir/${OutputSuffix}.${genome}.cufflinks.log
fi
if [ $runcufflinks == "1" ];then
  mv $CufflinksDir/genes.fpkm_tracking $CufflinksDir/${OutputSuffix}.${genome}.genes.fpkm_tracking
  mv $CufflinksDir/isoforms.fpkm_tracking $CufflinksDir/${OutputSuffix}.${genome}.isoforms.fpkm_tracking
  mv $CufflinksDir/transcripts.gtf $CufflinksDir/${OutputSuffix}.${genome}.transcripts.gtf
  echo "   Done running cufflinks"
else
  echo "   Skipping cufflinks..."
fi

echo "*****************************************************************************"
echo "8. Direct transcriptome mapping using Salmon:"
echo "   Salmon will automatically detect library type by default (-l A)"
if [ $mode == "single" ];then
echo "   Using strand type: "$StrandType
salmon quant \
    -i $HomeDir/$genome/Index/SalmonIndex \
    -p $CPU \
    -l $StrandType \
    -r $Data \
    -o $SalmonOutputDir \
                2> $SalmonOutputDir/${OutputSuffix}.${genome}.salmon.log
else
echo "   Using strand type: "$StrandType
salmon quant \
    -i $HomeDir/$genome/Index/SalmonIndex \
    -p $CPU \
    -l $StrandType \
    -1 $DataLeft \
    -2 $DataRight \
    -o $SalmonOutputDir \
                2> $SalmonOutputDir/${OutputSuffix}.${genome}.salmon.log
fi
echo "   Done Salmon quantification."
echo "   You can look at the strand info in lib_format_counts.json file"
mv $SalmonOutputDir/quant.sf $SalmonOutputDir/${OutputSuffix}.${genome}.quant.sf
mkdir $SalmonOutputDir/${OutputSuffix}.salmon
mv $SalmonOutputDir/aux_info $SalmonOutputDir/${OutputSuffix}.salmon/ && mv $SalmonOutputDir/cmd_info.json $SalmonOutputDir/${OutputSuffix}.salmon/
mv $SalmonOutputDir/lib_format_counts.json $SalmonOutputDir/${OutputSuffix}.salmon/ && mv $SalmonOutputDir/libParams $SalmonOutputDir/${OutputSuffix}.salmon/
mv $SalmonOutputDir/logs $SalmonOutputDir/${OutputSuffix}.salmon/

if [ -d $HomeDir/$genome/Index/SalmonIndexWithTE ];then
echo "8.1 Direct transcriptome mapping using Salmon, including TE:"
if [ $mode == "single" ];then
echo "   Using strand type: "$StrandType
salmon quant \
    -i $HomeDir/$genome/Index/SalmonIndexWithTE \
    -p $CPU \
    -l $StrandType \
    -r $Data \
    -o $SalmonOutputTEDir \
                2> $SalmonOutputTEDir/${OutputSuffix}.${genome}.salmonWithTE.log
else
echo "   Using strand type: "$StrandType
salmon quant \
    -i $HomeDir/$genome/Index/SalmonIndexWithTE \
    -p $CPU \
    -l $StrandType \
    -1 $DataLeft \
    -2 $DataRight \
    -o $SalmonOutputTEDir \
                2> $SalmonOutputTEDir/${OutputSuffix}.${genome}.salmonWithTE.log
fi
echo "   Done Salmon quantification, including TE."
mv $SalmonOutputTEDir/quant.sf $SalmonOutputTEDir/${OutputSuffix}.${genome}.WithTE.quant.sf
mkdir $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE
mv $SalmonOutputTEDir/aux_info $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE/ && mv $SalmonOutputTEDir/cmd_info.json $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE/
mv $SalmonOutputTEDir/lib_format_counts.json $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE/ && mv $SalmonOutputTEDir/libParams $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE/
mv $SalmonOutputTEDir/logs $SalmonOutputTEDir/${OutputSuffix}.salmonWithTE/
fi

echo "*****************************************************************************"
echo "9. Quantify gene abundance using featureCounts from Subread package:"
if [ $mode == "single" ];then
echo "   Running featureCounts, single-end mode"
echo "   featureCounts does not count reads overlapping with more than one feature"
echo "   feature_counts/*gene.txt file contains tidy gene counts"
featureCounts \
    -T $CPU \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.log

    if [ $pairrpf == "1" ];then
    featureCounts \
    -T $CPU \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.log
    
    featureCounts \
    -T $CPU \
    -t CDS \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.log
    fi
    
else
echo "   Running featurecounts, paired-end mode"
echo "   featureCounts does not count reads overlapping with more than one feature"
echo "   feature_counts/*.gene.txt file contains tidy gene counts"
featureCounts \
    -T $CPU \
    -p \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -D 2000 \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.log
    
    if [ $pairrpf == "1" ];then
    featureCounts \
    -T $CPU \
    -p \
    -t exon \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -D 2000 \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.log
    
    featureCounts \
    -T $CPU \
    -p \
    -t CDS \
    -g gene_id \
    -a $HomeDir/$genome/Annotation/${genome}.RefSeq.reduced.bed12.geneid.gtf \
    -D 2000 \
    -o $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.FullTable.txt \
    $GenomeMappingDir/${OutputSuffix}.${genome}.sorted.unique.bam \
    2> $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.log
    fi
fi
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.gene.txt
if [ $pairrpf == "1" ];then
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.FullTable.txt | awk 'BEGIN{OFS="\t"}{print $1,$6,$7}' > \
    $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt
awk '{print $1}' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.namesall
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt|awk '{print $1}' > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names
awk 'NR>1' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt|awk '{print $1}' > $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt.names
cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt.names|sort|uniq -c|awk '$1==1'|awk '{print $2}' > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names.lncRNA
for name in `cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names.lncRNA`;do awk -v t=$name '{if ($1==t) print $0}' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt;done > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.lncRNA
	
cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.lncRNA > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.All.CDS.gene.txt.unorder
for name in `cat $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.namesall`;do awk -v t=$name '{if ($1==t) print $0}' $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.All.CDS.gene.txt.unorder;done > \
	$FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.All.CDS.gene.txt
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.namesall $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.names.lncRNA
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.mRNA.CDS.gene.txt.names
mv $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.gene.txt.lncRNA $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.lncRNA.gene.txt
rm -rf $FeatureDir/${OutputSuffix}.${genome}.featureCounts.unique.All.CDS.gene.txt.unorder
fi

echo "   Done running featureCounts"

echo "*****************************************************************************"
echo "10. Finishing:"
if [ $bigWig == "1" ];then
Norm=`calculating $AllMapReads/1000000`
ScalingPlus=`calculating 1/$Norm`
ScalingMinus=`calculating -1/$Norm`
echo "   Getting normalization: "$Norm
Databam=$GenomeMappingDir/${OutputSuffix}.${genome}.sorted.bam
DatabamPre=$GenomeMappingDir/${OutputSuffix}.${genome}.sorted

if [ $mode == "single" ];then
echo "   Generating track files for single end data"
bedtools genomecov -bga -split -strand + -ibam ${Databam} -scale $ScalingPlus > ${DatabamPre}.plus.bedGraph
bedtools genomecov -bga -split -strand - -ibam ${Databam} -scale $ScalingMinus > ${DatabamPre}.minus.bedGraph
bedtools genomecov -bga -split -ibam ${Databam} -scale $ScalingPlus > ${DatabamPre}.bedGraph
awk '$4!=0' ${DatabamPre}.plus.bedGraph > ${DatabamPre}.plus.filtered.bedGraph
awk '$4!=0' ${DatabamPre}.minus.bedGraph > ${DatabamPre}.minus.filtered.bedGraph
awk '$4!=0' ${DatabamPre}.bedGraph > ${DatabamPre}.filtered.bedGraph
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.plus.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.plus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.minus.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.minus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.bedGraph.bw
rm -rf ${DatabamPre}.plus.bedGraph ${DatabamPre}.minus.bedGraph ${DatabamPre}.bedGraph
rm -rf ${DatabamPre}.plus.filtered.bedGraph ${DatabamPre}.minus.filtered.bedGraph ${DatabamPre}.filtered.bedGraph
else
echo "   Generating track files for paired end data"
samtools view -b ${Databam}|bedtools bamtobed -bed12 > ${DatabamPre}.bed12
awk 'BEGIN{OFS="\t"}{strand=$6;if (strand=="+") revstrand="-";else if (strand=="-") revstrand="+";if ($4~/1$/) {$6=revstrand;print $0;} else print $0}' ${DatabamPre}.bed12 > ${DatabamPre}.bed12.rev
bedtools genomecov -bga -split -strand + -i ${DatabamPre}.bed12.rev -scale $ScalingPlus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DatabamPre}.plus.bedGraph
bedtools genomecov -bga -split -strand - -i ${DatabamPre}.bed12.rev -scale $ScalingMinus -g $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt > ${DatabamPre}.minus.bedGraph
bedtools genomecov -bga -split -ibam ${Databam} -scale $ScalingPlus > ${DatabamPre}.bedGraph
awk '$4!=0' ${DatabamPre}.plus.bedGraph > ${DatabamPre}.plus.filtered.bedGraph
awk '$4!=0' ${DatabamPre}.minus.bedGraph > ${DatabamPre}.minus.filtered.bedGraph
awk '$4!=0' ${DatabamPre}.bedGraph > ${DatabamPre}.filtered.bedGraph
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.plus.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.plus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.minus.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.minus.bedGraph.bw
$HomeDir/bin/bedGraphToBigWig ${DatabamPre}.filtered.bedGraph $HomeDir/$genome/Sequence/${genome}.ChromInfo.txt ${DatabamPre}.bedGraph.bw
rm -rf ${DatabamPre}.plus.bedGraph ${DatabamPre}.minus.bedGraph ${DatabamPre}.bedGraph
rm -rf ${DatabamPre}.plus.filtered.bedGraph ${DatabamPre}.minus.filtered.bedGraph ${DatabamPre}.filtered.bedGraph
rm -rf ${DatabamPre}.bed12.rev ${DatabamPre}.bed12
fi

mv ${DatabamPre}*bedGraph.bw $TracksDir

fi

echo "Time Started: "$St
Ed=`date`
echo "Time Ended:   "$Ed
echo "*                           End of the pipeline                             *"
echo "*****************************************************************************"
