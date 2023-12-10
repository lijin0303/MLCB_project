#!/bin/bash
#SBATCH -n 11
#SBATCH -N 1
#SBATCH -p sched_any,newnodes
#SBATCH -t 0-12:00:00
#SBATCH -o calling_%j.o
#SBATCH -e calling_%j.o
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rachel95912@gmail.com
tc=10

workdir=/home/ruitongl/MLCB_project/peak
source /home/ruitongl/MLCB_project/phasecall/bin/activate

python ${workdir}/cPeaks/main.py --fragment_path ${workdir}/fragments.sorted.tsv.gz \
      --barcode_path ${workdir}/atac_CB.txt --output ${workdir}/SNU601_peak \
      --num_cores ${tc}