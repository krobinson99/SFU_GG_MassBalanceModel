#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --account=def-gflowers
#SBATCH --mem=3G
#SBATCH -o calculate_mb-%A_%a.out

#SBATCH --array=0

echo calculate glacier-wide net mass balance
python calculate_mb.py $SLURM_ARRAY_TASK_ID
echo done!
