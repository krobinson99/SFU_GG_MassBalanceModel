#!/usr/bin/env bash

module load cdo

# Raw NARR files are kept in scratch to save space
cd /home/krobin/scratch/RawNARR/TEST

# merge all monthly air and hgt files into 1 file per year (apcp is already in the right format)
for year in {1979..2022}
do
cdo mergetime air.${year}{01..12}.nc air.${year}.nc
cdo mergetime hgt.${year}{01..12}.nc hgt.${year}.nc
done

# next cut files down to KW domain (saves storage space)
for year in {1979..2022}
do
cdo sellonlatbox,-140.5,-138.,60.,61.5 air.${year}.nc KRH_air.${year}.nc
cdo sellonlatbox,-140.5,-138.,60.,61.5 hgt.${year}.nc KRH_hgt.${year}.nc
cdo sellonlatbox,-140.5,-138.,60.,61.5 apcp.${year}.nc KRH_pcp.${year}.nc
done

# copy the NARR inputs to some folder (ideally not scratch) where they can be used as input to the downscaling scripts
cp KRH_air.* /home/krobin/projects/def-gflowers/krobin/NARRinputs
cp KRH_hgt.* /home/krobin/projects/def-gflowers/krobin/NARRinputs
cp KRH_apcp.* /home/krobin/projects/def-gflowers/krobin/NARRinputs


