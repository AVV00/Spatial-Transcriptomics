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
#For PE data, the adapter sequence auto-detection is disabled by default since the adapters can be trimmed by overlap analysis.
# However, you can specify --detect_adapter_for_pe to enable it.

fastp_cmd="$fastp -i $inputDir/$sampleId*R1*.gz -I $inputDir/$sampleId*R2*.gz \
-o $cleanFastq/$sampleId.R1.clean.fq.gz -O $cleanFastq/$sampleId.R2.clean.fq.gz \
-j $cleanFastq/$sampleId.json -h $cleanFastq/$sampleId.html \
--qualified_quality_phred 20 --unqualified_percent_limit 30 \
--n_base_limit 5 --detect_adapter_for_pe;"

echo $fastp_cmd

echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ fastp running ------\n"

if [ -z $4 ] || [ $4 = 0 ];then
   eval $fastp_cmd
else
   echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ skip fastp running ------\n"
fi

echo -e `date +"%Y-%m-%d %H:%M:%S"` "\n------ fastp finished ------\n"


