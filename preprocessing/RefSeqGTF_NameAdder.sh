#!/usr/bin/env bash
#RefSeqGTF_NameAdder.sh
#This is a pipeline to substitue column 10 to be the read gene names.
#Required input: Any gtf file with column 10 = RefSeq ID, And a matching file from current values to new values.
#Example of the matching file: NM_001000000 Olr1389
#Version: Yu Sun, 2022-01-28

if [ ! -n "$1" ]
then
  echo "    This is a pipeline to substitue column 10 to be the read gene names."
  echo "    Required input: Any gtf file with column 10 = RefSeq ID, And a matching file from current values to new values."
  echo "    Usage: `basename $0` [Ori.gtf] [Unique Matching|oriName newName] [New.gtf]"
  echo "    Output: a new gtf file"
else
  Data=$1
  Matching=$2
  NewFile=$3
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' $Data > $Data.1-9
  awk '{print $11,$12,$13,$14,$15,$16}' $Data> $Data.11-16
  awk '{print $10}' $Data |sed 's/"//g'|sed 's/;//g' > $Data.transcript
  for name in `cat $Data.transcript`;do awk -v getnames=$name '{if ($1==getnames) print "\""$2"\""";"}' $Matching;done > $Data.gene
  paste $Data.1-9 $Data.gene $Data.11-16 -d" " > $NewFile
  rm -rf $Data.1-9 && rm -rf $Data.11-16 && rm -rf $Data.transcript && rm -rf $Data.gene
fi
