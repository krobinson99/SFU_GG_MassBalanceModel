#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=def-gflowers
#SBATCH --mem=40G
#SBATCH -o mbm-%A_%a.out

#SBATCH --array=0

echo starting mass-balance model # saving all variables: 4 hr, 40G
python massbalance.py $SLURM_ARRAY_TASK_ID
echo done!
