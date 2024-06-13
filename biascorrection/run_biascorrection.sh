#!/bin/bash
#SBATCH --time=0:10:00
#SBATCH --account=def-gflowers
#SBATCH --mem=15G
#SBATCH -o biascorrect-%A_%a.out

#SBATCH --array=2001


echo bias correcting NARR temperature and precipitation
python biascorrection.py $SLURM_ARRAY_TASK_ID
echo done!
