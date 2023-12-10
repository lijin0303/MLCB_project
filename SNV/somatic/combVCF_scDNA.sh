#!/bin/bash

path="/home/ruitongl/MLCB_project/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

bcftools=Monopogen/apps/bcftools
bgzip=Monopogen/apps/bgzip
cd /home/ruitongl/MLCB_project

vcft=gl
c=1
while(($c < 19))
do
    echo $c
    find DNA_germline -type f -name chr${c}:*.${vcft}.vcf.gz -exec ${bcftools} index {} -f -t 4  \;
    ${bcftools} concat DNA_germline/chr${c}:*.${vcft}.vcf.gz -Ou | ${bcftools} sort -Ov -o DNA_germlineC/chr${c}.${vcft}.vcf
    ${bgzip} -c DNA_germlineC/chr${c}.${vcft}.vcf > DNA_germlineC/chr${c}.${vcft}.vcf.gz # -Oz does not properly compress the file
    rm DNA_germlineC/chr${c}.${vcft}.vcf
    ${bcftools} index DNA_germlineC/chr${c}.${vcft}.vcf.gz -f -t 4
    ((c=c+1))
done

vcft=phased
c=1
while(($c < 19))
do
    echo $c
    find DNA_germline -type f -name chr${c}:*.${vcft}.vcf.gz -exec ${bcftools} index {} -f -t 4  \;
    ${bcftools} concat DNA_germline/chr${c}:*.${vcft}.vcf.gz -Ou | ${bcftools} sort -Ov -o DNA_germlineC/chr${c}.${vcft}.vcf
    ${bgzip} -c DNA_germlineC/chr${c}.${vcft}.vcf > DNA_germlineC/chr${c}.${vcft}.vcf.gz # -Oz does not properly compress the file
    rm DNA_germlineC/chr${c}.${vcft}.vcf
    ${bcftools} index DNA_germlineC/chr${c}.${vcft}.vcf.gz -f -t 4
    ((c=c+1))
done



# cp /pool001/meganle/DNAseq_results/germline/*.gl.vcf.gz /home/ruitongl/MLCB_project/DNA_germline/