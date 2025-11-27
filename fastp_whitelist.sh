#!/bin/bash

#fastp.sh
#this script is for running fastq for qc and trimming on fastq files
#usage 
#--- 
# ./fastp.sh sampleId inputDir outputDir 

###################################################################

# Usage help and parameters
if [ $# -lt 3 ];then
  echo "Usage: $0 <sampleId> <inputDir:absPath> <outputDir:absPath> [skip:0 | 1]"
  exit 1
fi

sampleId=$1
inputDir=$2
outputDir=$3

# softwares and dir
fastp=/public/home/miaozhu_gibh/software/fastp
cleanFastq=$outputDir/cleanFastq


# Prepare: change to workspace
if [ ! -d $cleanFastq ]; then
  mkdir -p $cleanFastq
fi

cd $cleanFastq

#prepare command

fastp_cmd="$fastp -i $inputDir/$sampleId*.gz  \
-o $cleanFastq/$sampleId.clean.fq.gz \
--qualified_quality_phred 20 --unqualified_percent_limit 30 \
--n_base_limit 1 --length_required 30 \
-A -G \
-j $cleanFastq/$sampleId.json -h $cleanFastq/$sampleId.html ;"

echo $fastp_cmd

echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ fastp running ------\n"

if [ -z $4 ] || [ $4 = 0 ];then
   eval $fastp_cmd
else
   echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ skip fastp running ------\n"
fi

echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ fastp finished ------\n"


