#!/bin/bash
i=1
while(($i < 23))
do
    sbatch submit_get_imputation $i
    ((i=i+1))
done

sbatch submit_get_imputation X
