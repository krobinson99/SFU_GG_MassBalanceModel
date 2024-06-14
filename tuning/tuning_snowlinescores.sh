#!/bin/bash
#SBATCH --time=00:40:00
#SBATCH --account=def-gflowers
#SBATCH --mem=30G
#SBATCH -o calculate_ss-%A_%a.out

#SBATCH --array=0

echo calculate snowline scores for each delineated satellite image
python calculate_snowlinescore.py  $SLURM_ARRAY_TASK_ID
echo done!
