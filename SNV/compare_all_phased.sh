#!/bin/sh

#Report the comparion of germline SNVs between DNAseq and ATACseq

for FILE in /pool001/meganle/scATAC_results/germline/*.phased.vcf.gz
do 
name=${FILE##*/} #"Shell parameter expansion"
region=${name%%.*}
chr=${region%:*}
#echo $FILE $name $region $chr
vcftools --gzvcf  $FILE  --gzdiff /pool001/meganle/DNAseq_results/germline/${name} --diff-discordance-matrix --out /pool001/dschaffe/variant-calling/germline_comparison/${region}_phased-phased --chr $chr
done
