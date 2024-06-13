#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --account=def-gflowers
#SBATCH --mem=2G

#download big NARR files in scratch directory - they don't need to be kept forever
cd /home/krobin/scratch/RawNARR/TEST

# download precip data (one file per year with daily prcp values):
for year in {1979..2022}
do
wget -nc -c -nv "https://psl.noaa.gov/thredds/fileServer/Datasets/NARR/Dailies/monolevel/apcp.${year}.nc"
done

#download air temperature data and geopotential heigh data (one file per month with 3-hourly values):
for year in {1979..2022}
do
for month in {01..12}
do
wget -nc -c -nv "https://downloads.psl.noaa.gov/Datasets/NARR/pressure/air.${year}${month}.nc"
wget -nc -c -nv "https://downloads.psl.noaa.gov/Datasets/NARR/pressure/hgt.${year}${month}.nc"
done
done


