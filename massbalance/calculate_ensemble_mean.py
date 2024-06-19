# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:32:38 2023

This script takes N simulations that comprise the model ensemble (generally N=100) 
and calculates the ensemble mean and standard deviation for each year and each model output variable. 

@author: katierobinson
"""

import numpy as np
import os
import sys
import pandas as pd
from netCDF4 import Dataset
# Import parameters from config file
from mbm_config import Easting_grid, Northing_grid, Model_functions
# Import functions:
sys.path.insert(1,Model_functions)
from model_functions import save_to_netcdf

# =============================================================================
# Paths to model outputs:
# =============================================================================

Glacier_ID = 'KRH'                                                               # identifier for the glacier. All inputs/outputs should also have the same ID.
model_runs = '/home/krobin/scratch/ModelRuns/KRH/VAAtest1'                       # Directory where model outputs are stored 
output_folder = '/home/krobin/scratch/ModelRuns/KRH/VAAtest1/ConcatenatedRuns'   # Directory where the ensemble means will be saved 
N_sims = 100                                                                     # number of simulations in the model ensemble 
OUTPUT_ID = 99999                                                                # use an integer > N_sims to distinguish the ensemble mean from the individual runs

# =============================================================================
# Get year from the job array
# =============================================================================

sims = np.arange(0,N_sims) 

year = int(sys.argv[1])
print('cancatenating runs for ',year)
print(type(year))

# =============================================================================
# Load input geometry:
# =============================================================================

Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)

# =============================================================================
# Function to calculate the ensemble mean and standard deviation
# =============================================================================
def concatenate_model_ouputs(year, sims, varname, var):
    '''
    Inputs: 
        year
        simulations (array or list of integers corresponding to the model runs to be averaged: e.g. 1--100)
        varname (variable name string, e.g. 'Icemelt')
        var: netcdf variable (string, e.g. 'Ice melt')
    Returns:
        Mean and standard deviation of all simulations for a given year and variable. 
    '''
    
    print('year:',year)
    sys.stdout.flush()
    dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')  
    
    # =============================================================================
    #   Calculate mean of all simulation ensemble:
    # =============================================================================
    all_sims = np.empty((len(sims),len(dates),Xgrid.shape[0],Xgrid.shape[1]))
    for sim in sims:
        print('Calculating mean and std dev. for var:',str(var),str(sim))
        sys.stdout.flush()
        
        inMB = Dataset(os.path.join(model_runs,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
        mb = inMB.variables[var][:]
        sys.stdout.flush()
        inMB.close()
        
        # Create list of all data for calculating the std dev.
        all_sims[sim,:,:,:] = mb
        
    mean_mb = np.nanmean(all_sims,axis=0)
    std = np.std(all_sims,axis=0)
        
    save_to_netcdf(mean_mb, var, os.path.join(output_folder,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(OUTPUT_ID) + '.nc'), year, Xgrid, Ygrid) 
    save_to_netcdf(std, var, os.path.join(output_folder,varname + '_' + 'std' + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(OUTPUT_ID) + '.nc'), year, Xgrid, Ygrid) 
    
    
# Call function for each of the model output variables:    
# =============================================================================
concatenate_model_ouputs(year,sims,'Snowmelt', 'Snowmelt',24)   
concatenate_model_ouputs(year,sims,'Refrozenmelt', 'Refrozenmelt',24)   
concatenate_model_ouputs(year,sims,'Netsnowmelt', 'Netsnowmelt',24)  
concatenate_model_ouputs(year,sims,'Glaciericemelt', 'Glaciericemelt',24)   
concatenate_model_ouputs(year,sims,'Superimposedicemelt', 'Superimposedicemelt',24)  
concatenate_model_ouputs(year,sims,'Rain', 'Rain',24)    
concatenate_model_ouputs(year,sims,'Refrozenrain', 'Refrozenrain',24)  
concatenate_model_ouputs(year,sims,'Rainrunoff', 'Rainrunoff',24) 
concatenate_model_ouputs(year,sims,'Accumulation', 'Accumulation',24)
concatenate_model_ouputs(year,sims,'Snowdepth', 'Snowdepth',24) 
concatenate_model_ouputs(year,sims,'Ptau', 'Ptau',24)
concatenate_model_ouputs(year,sims,'Superimposedice', 'Superimposedice',24)  
concatenate_model_ouputs(year,sims,'Netbalance', 'Netbalance',24)
# =============================================================================

  
