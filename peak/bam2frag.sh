#!/bin/bash
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -p sched_any,newnodes
#SBATCH -t 0-12:00:00
#SBATCH -o calling_%j.o
#SBATCH -e calling_%j.o
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rachel95912@gmail.com
tc=5

workdir=/home/ruitongl/MLCB_project/peak
source /home/ruitongl/MLCB_project/phasecall/bin/activate
bgzip=/home/ruitongl/MLCB_project/Monopogen/apps/bgzip
tabix=/home/ruitongl/MLCB_project/Monopogen/apps/tabix

# ATACbam=/pool001/meganle/SNU601_scATAC.bam 
path="/home/ruitongl/MLCB_project/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

sinto fragments -b ${ATACbam} -p ${tc} \
     -f ${workdir}/fragments.tsv \
     -c ${workdir}/atac_CB.txt

# sort, compress, and index
sort -k1,1 -k2,2n ${workdir}/fragments.tsv > ${workdir}/fragments.sorted.tsv
${bgzip} -@ ${tc} ${workdir}/fragments.sorted.tsv
${tabix} -p bed ${workdir}/fragments.sorted.tsv.gz


