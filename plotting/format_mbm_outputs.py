# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:12:34 2024

This script formats mass-balance model outputs (shape = 365 x 230 x 329) as either 
a timeseries, averaged over the catchment or glacier area (shape = 365 x 1) or as 
an annual distributed field (shape = 1 x 230 x 329). 

Saving the formatted outputs as text files reduces the memory needed to store 
outputs from about 105 MB to 1MB per file.

This script also makes figure generation more efficient as the outputs are already
saved in the correct format needed for plotting.

@author: katierobinson
"""

import numpy as np

from plotting_functions import process_model_timeseries, save_timeseries, \
process_model_distributed_fields, save_distributed_fields

# =============================================================================
# paths and params for post-processing mass-balance model outputs:
# =============================================================================

model_outputs = 'D:/Model Runs/REF_MODEL/Sim_99999_v2'     # directory with the model outputs
Glacier_ID = 'KRH'                                         # identifier for the glacier. All inputs/outputs will also have the same ID.
sim = 99999                                                # sim ID of the model outputs that are being processed
processed_outputs = 'D:/Model Runs/REF_MODEL/Processing' # directory where processed outputs will be saved

# time period:
start_year = 1979
end_year = 1980

# Input geometry
Sfc_grid_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_SfcType.txt'    # text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.
KW_grid_file = 'F:/Mass Balance Model/Kaskawulsh-Mass-Balance/RunModel/KRH_Tributaries.txt' # text file where Kaskawulsh Glacier is finite, otherwise NaN 

# =============================================================================
# Set up glacier grids and variables for post-processing
# =============================================================================

Sfc = np.loadtxt(Sfc_grid_file)         # array where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.
KW_grid = np.loadtxt(KW_grid_file)      # array where non-NaN = Kaskawulsh Glacier, otherwise NaN 

glacierized_area_grid = np.array(Sfc)   # array where 0 = glacierized area, otherwise NaN
glacierized_area_grid[np.where(Sfc==1)] = np.nan

# variables names that will be used to save text files generated in this script
varnames = ['massbal','totalsnowmelt','refreezing','netsnowmelt','glaciericemelt','superimposedice_melt','rain','refrozen_rain','rain_runoff','accumulation','snowdepth','Ptau','SI']

# years to be processed:
years = np.arange(start_year,end_year)

# =============================================================================
# Calculate the catchment-wide average annual mass balance timeseries. Returns a daily timeseries for each hydrological year (Oct 1-Sep 30) in the array 'years', 
# =============================================================================

# Calculate the daily catchment-wide average of each variable:
catchment_wide_avg_timeseries = process_model_timeseries(model_outputs, Glacier_ID, sim, years, Sfc, False)

# Save the catchment-wide average timeseries for each variable as text files:
save_timeseries(list(catchment_wide_avg_timeseries),years,processed_outputs,varnames,'REF_MODEL','krh')

# Calculate the catchment-wide average standard deviation from the ensemble mean:
catchment_wide_avg_timeseries_stddev = process_model_timeseries(model_outputs, Glacier_ID, sim, years, Sfc, True)

# Save the standard deviation timeseries as text files:
save_timeseries(list(catchment_wide_avg_timeseries_stddev),years,processed_outputs,varnames,'REF_MODEL','krh_std')

# =============================================================================
# Repeat, but averaged over the Kaskawulsh Glacier only (not catchment-wide) 
# =============================================================================

# Calculate the daily glacier-wide average of each variable:
# Note: Here, Sfc is replaced with KW_grid, which defines the new area to be averaged over.
glacier_wide_avg_timeseries = process_model_timeseries(model_outputs, Glacier_ID, sim, years, KW_grid, False)

# Save the glacier-wide average timeseries for each variable as text files:
save_timeseries(list(glacier_wide_avg_timeseries),years,processed_outputs,varnames,'REF_MODEL','kw')

# Calculate the glacier-wide average standard deviation from the ensemble mean:
glacier_wide_avg_timeseries_stddev = process_model_timeseries(model_outputs, Glacier_ID, sim, years, KW_grid, True)

# Save the standard deviation timeseries as text files:
save_timeseries(list(glacier_wide_avg_timeseries_stddev),years,processed_outputs,varnames,'REF_MODEL','kw_std')

# =============================================================================
# Calculate the annual distributed fields for each variable
# =============================================================================

# compute and save annual distributed mass balance fields:
distributed_fields = process_model_distributed_fields(model_outputs, Glacier_ID, sim, years, False)
save_distributed_fields(list(distributed_fields),varnames,years,processed_outputs,'REF_MODEL',False)

# compute and save annual distributed standard deviation fields:
distributed_fields_std = process_model_distributed_fields(model_outputs, Glacier_ID, sim, years, True)
save_distributed_fields(list(distributed_fields_std),varnames,years,processed_outputs,'REF_MODEL',True)







