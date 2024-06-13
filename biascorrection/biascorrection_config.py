# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

#domain
Glacier_ID = 'KRH'                                          # string, used to name output files
Easting_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Xgrid.txt'   # path to text file containing array of easting coords for downscaled grid
Northing_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Ygrid.txt'  # path to text file containing array of northing coords for downscaled grid
Elev_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/Zgrids'           # path to directory containing text file with elevations for downscaled grid

# physics
R2S = 1.0 # rain to snow threshold (degree C)

# inputs
Model_functions = '/home/krobin/projects/def-gflowers/krobin/MBM/model_functions'                                         # path to directory with model functions
downscaled_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/DownscaledClimateInputs'       # path to directory with coarse NARR data (temp, precip, geopotential height)

# outputs
output_path = '/home/krobin/scratch/MBM_testing/BC_test' # path to save downscaled outputs
