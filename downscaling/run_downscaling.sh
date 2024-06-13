#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=def-gflowers
#SBATCH --mem=8G
#SBATCH -o downscale-%A_%a.out

#SBATCH --array=1999-2000

echo downscaling NARR temperature and precipitation
python downscaling.py $SLURM_ARRAY_TASK_ID
echo done!
