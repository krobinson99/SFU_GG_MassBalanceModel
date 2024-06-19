#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --account=def-gflowers
#SBATCH --mem=80G
#SBATCH -o Concat-%A_%a.out

#SBATCH --array=1979-2022

echo calculating ensemble mean and standard deviations
python calculate_ensemble_mean.py $SLURM_ARRAY_TASK_ID # for saving all vars at once use 80GB, 5hrs
echo done!
