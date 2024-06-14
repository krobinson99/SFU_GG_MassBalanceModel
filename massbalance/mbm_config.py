# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:15:08 2021

@author: katierobinson
"""

# Time control
start_year = 1979
end_year = 2022
timestep = 3           # timestep of input data (hours).

# Glacier id
Glacier_ID = 'KRH'     # identifier for the glacier. All inputs/outputs will also have the same ID.

# Melt model parameters
params_file = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/Params/REF_MODEL/Tuning_Params_Final_REFMODEL.csv' # path to csv file containing aice, asnow, MF

# physics
R2S = 1.0              # rain-to-snow threshold (degree C).

# Input data
Model_functions = '/home/krobin/projects/def-gflowers/krobin/MBM/model_functions'                                   # path to directory with model functions
Precip_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/BiasCorrectedClimateInputs'  # path to directory with precip inputs
Temp_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/BiasCorrectedClimateInputs'    # path to directory with temperature inputs
Solar_inputs = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/SolarInputs'                  # path to directory with PDCSR inputs

# Input geometry
Easting_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Xgrid.txt'   # array of easting coords for downscaled grid
Northing_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_Ygrid.txt'  # array of northing coords for downscaled grid
Sfc_grid = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Downscaling/InputGeometry_KRH/KRH_SfcType.txt'     # array that determines the surface type of each gridcell in domain (1 = off-glacier, 0 = on-glacier, NaN = not in the domain).

# Debris parameters
debris_parameterization = 'Sub-debris melt scaling'                                                                    # Debris treatment. Options are (1) None (no debris), (2) 'Boolean debris' (sets a_ice = 0 in debris cells), or (3) 'Sub-debris melt scaling' (thickness dependent scaling)
debris_map = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/InputGeometries/KRH_debrismap.txt' # Path to text file. If debris parameterization is (1) file should be same as Sfc_grid; (2) debris cells are finite, else NaN, (3) cells w debris have value = debris thickness, else NaN
cleanice_melt = 2.0277           # clean ice melt (m w.e).
peak_melt = 2.1717               # peak melt (m w.e.).
peak_melt_thickness = 0.006      # thickness at which peak melt occurs (m).
transition_thickness = 0.019     # thickness at which melt equals clean ice melt (m) (aka critical thickness).
b0 = 11.0349260206858            # b0 and k are site-specific params from Rounce et al. (2021).
k = 1.98717418666925

# Outputs
OUTPUT_PATH = '/home/krobin/scratch/ModelRuns/KRH/test_2024_06_14'  # directory where model outputs should be saved
output_runoff = True         # saves glacier ice melt, net snowmelt, superimposed ice melt, runoff from rain
output_refreezing = True     # saves total snowmelt, refrozen snowmelt, total rainfall, refrozen rainfall, Ptau, and the superimposed ice layer
output_accumulation = True   # saves accumulation and the snowpack layer.
output_massbalance = True    # saves mass balance
