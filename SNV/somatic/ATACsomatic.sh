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
python=/home/ruitongl/MLCB_project/phasecall/bin/python
path="/home/ruitongl/MLCB_project/Monopogen"
CBcnt="/home/ruitongl/MLCB_project/readcnt/SNU601_scATAC_scReadCnt.csv"
region="/home/ruitongl/MLCB_project/refdata/GRCh38.chr$1.region.lst"
fa="/pool001/meganle/hs38.fa"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

module add R/3.6.2

fileP=/pool001/ruitongl/scATAC_results_$1/

mkdir /pool001/ruitongl/scATAC_results_$1/
mkdir /pool001/ruitongl/scATAC_results_$1/Bam/
mkdir /pool001/ruitongl/scATAC_results_$1/germline/
mkdir /pool001/ruitongl/scATAC_results_$1/Script/

cp -r /pool001/meganle/scATAC_results/Bam/*chr$1.* /pool001/ruitongl/scATAC_results_$1/Bam/
cp 	germlineC/chr$1.* 	/pool001/ruitongl/scATAC_results_$1/germline/


${python}  ${path}/src/Monopogen2.py  somatic \
    -a   ${path}/apps  -r  ${region}  -t ${tc} \
	-i  ${fileP}  -l  ${CBcnt}   -s featureInfo \
	-g  ${fa}

${python}  ${path}/src/Monopogen2.py  somatic \
	-a   ${path}/apps  -r  ${region}  -t ${tc}  -w 10MB \
	-i  ${fileP}  -l  ${CBcnt} -s cellScan \
	-g   ${fa}

${python}  ${path}/src/Monopogen2.py  somatic  \
	-a   ${path}/apps  -r  ${region}  -t ${tc} \
	-i  ${fileP}  -l  ${CBcnt}   -s LDrefinement \
	-g   ${fa}


cp /pool001/ruitongl/scATAC_results_$1/somatic/chr$1.cell.txt /home/ruitongl/SNV_final/CB/
cp /pool001/ruitongl/scATAC_results_$1/somatic/chr$1.SNV_mat.RDS /home/ruitongl/SNV_final/mutmat/
cp /pool001/ruitongl/scATAC_results_$1/somatic/chr$1.putativeSNVs.csv /home/ruitongl/SNV_final/SNVinfo/
rm -r /pool001/ruitongl/scATAC_results_$1/
