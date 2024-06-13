# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

# time control
time_step = 3                                               # timestep of NARR temperature data (hours)

#domain
Glacier_ID = 'KRH'                                          # string, used to name output files
UTM = 7                                                     # UTM zone
NARR_subregions = [9,10,14,15,16,19,20,21,22,25,26,27,28]   # indices of subregions used for precip downscaling, picked manually to omit points on opposite side of divide

# inputs
Model_functions = '/home/krobin/projects/def-gflowers/krobin/MBM/model_functions'                                         # path to directory with model functions
Climate_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/CoarseNARR_KRH'                  # path to directory with coarse NARR data (temp, precip, geopotential height)
Coarse_DEM_input = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/CoarseNARR_KRH/KRH_CE.nc'      # path to netcdf file with elevations of coarse NARR gridcells. From Datasets/NARR/time_invariant/hgt.sfc.nc (https://psl.noaa.gov/thredds/catalog/Datasets/NARR/time_invariant/catalog.html?dataset=Datasets/NARR/time_invariant/hgt.sfc.nc)
Easting_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Xgrid.txt'   # path to text file containing array of easting coords for downscaled grid
Northing_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Xgrid.txt'  # path to text file containing array of northing coords for downscaled grid
Elev_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/Zgrids'           # path to directory containing text file with elevations for downscaled grid

# outputs
OUTPUT_PATH = '/home/krobin/scratch/MBM_testing/Downscaling_test' # path to save downscaled outputs
