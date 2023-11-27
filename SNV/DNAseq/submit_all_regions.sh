#!/bin/bash

i=1

while(($i < 12))
do
    echo $i
    sbatch submit_calling $i
    ((i=i+1))
done

echo "12-13"
sbatch submit_calling 12-13
echo "14-15"
sbatch submit_calling 14-15
echo "16-18"
sbatch submit_calling 16-18
echo "19-22"
sbatch submit_calling 19-22
