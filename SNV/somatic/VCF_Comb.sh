#!/bin/bash
#SBATCH -n 5
#SBATCH -N 1
#SBATCH -p sched_any
#SBATCH -t 0-02:00:00
#SBATCH -o calling_%j.o
#SBATCH -e calling_%j.o
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rachel95912@gmail.com

bcftools=Monopogen/apps/bcftools
bgzip=Monopogen/apps/bgzip
cd /home/ruitongl/MLCB_project

vcft=gl
c=1
while(($c < 19))
do
    echo $c
    find germline -type f -name chr${c}:*.${vcft}.vcf.gz -exec ${bcftools} index {} -f -t 4  \;
    ${bcftools} concat germline/chr${c}:*.${vcft}.vcf.gz -Ou | ${bcftools} sort -Ov -o germlineC/chr${c}.${vcft}.vcf
    ${bgzip} -c germlineC/chr${c}.${vcft}.vcf > germlineC/chr${c}.${vcft}.vcf.gz # -Oz does not properly compress the file
    rm germlineC/chr${c}.${vcft}.vcf
    ${bcftools} index germlineC/chr${c}.${vcft}.vcf.gz -f -t 4
    ((c=c+1))
done

vcft=phased
c=1
while(($c < 19))
do
    echo $c
    find germline -type f -name chr${c}:*.${vcft}.vcf.gz -exec ${bcftools} index {} -f -t 4  \;
    ${bcftools} concat germline/chr${c}:*.${vcft}.vcf.gz -Ou | ${bcftools} sort -Ov -o germlineC/chr${c}.${vcft}.vcf
    ${bgzip} -c germlineC/chr${c}.${vcft}.vcf > germlineC/chr${c}.${vcft}.vcf.gz # -Oz does not properly compress the file
    rm germlineC/chr${c}.${vcft}.vcf
    ${bcftools} index germlineC/chr${c}.${vcft}.vcf.gz -f -t 4
    ((c=c+1))
done



# cp /pool001/meganle/DNAseq_results/germline/*.gl.vcf.gz /home/ruitongl/MLCB_project/DNA_germline/