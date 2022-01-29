#!/usr/bin/env bash
#RefSeqGTF_Simplifier.sh
#This is a pipeline to simply downloaded RefSeq.gtf file, remove _dup transcripts, and remove duplication of same transcript names on different chromosomes.
#Version: Yu Sun, 2022-01-28

if [ ! -n "$1" ]
then
  echo "    This is a pipeline to simply downloaded RefSeq.gtf file, remove _dup transcripts, and remove duplication of same transcript names on different chromosomes."
  echo "    Usage: `basename $0` [Genome.gtf]"
  echo "    Output: Genome.simplified.gtf"
else
  Data=$1
  rm -rf ${Data%.*}.simplified.gtf
  grep -v _dup $Data > ${Data%.*}.nodup.gtf
  awk 'BEGIN{OFS="\t"}{print $10,$1}' ${Data%.*}.nodup.gtf|sort -k1,1|uniq|sort -k1,1 -u > ${Data%.*}.nodup.chrmatchlist
  Lines=`wc -l ${Data%.*}.nodup.chrmatchlist|awk '{print $1}'`
  echo $Lines
  for i in $(seq 1 $Lines)
  do
    #echo $i
    Curr=`head -$i ${Data%.*}.nodup.chrmatchlist|tail -1`
    #echo $Curr
    Chr=`head -$i ${Data%.*}.nodup.chrmatchlist|tail -1|awk '{print $2}'`
    Trx=`head -$i ${Data%.*}.nodup.chrmatchlist|tail -1|awk '{print $1}'`
    awk -v Chrr=$Chr -v Trxx=$Trx '{if ($1==Chrr && $10==Trxx) print $0}' ${Data%.*}.nodup.gtf >> ${Data%.*}.simplified.gtf
  done
fi
