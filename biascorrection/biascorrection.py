# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 09:22:23 2021

Bias correct downscaled NARR temperature and precipitation

Inputs:
    Downscaled NARR T/P
Outputs:
    Bias corrected T/P
    
Notes:
- Temperature bias correction is the time-dependent factors from Young et al. (2021)
- Accumulation bias correction is elevation-dependent factors from Robinson et al. (2024) (see https://summit.sfu.ca/item/38185 for details)
- Only accumulation is bias corrected (precipitation that occurs when T > rain-to-snow threshold)

@author: katierobinson
"""
import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import sys
import os

# Import parameters from config file
from biascorrection_config import Glacier_ID, R2S
from biascorrection_config import downscaled_inputs, Easting_grid, Northing_grid, Elev_inputs, Model_functions, output_path

# Import functions for model:
sys.path.insert(1,Model_functions)
from model_functions import write_config_file, save_to_netcdf

# Save configuration file for this run to output directory:
write_config_file(output_path,"biascorrection_config.py")


# TEMPERATURE BIAS CORRECTION
# =============================================================================
DT = interp1d(np.asarray([1.,32.,60.,91.,121.,152.,182.,213.,244.,274.,305.,335.,366.]),np.asarray([-4.86,-3.98,-2.1,0.03,0.89,1.42,0.94,0.32,-0.32,-1.72,-4.37,-5.22, -4.86]),kind = 'linear')   
# =============================================================================  

# ACCUMULATION BIAS CORRECTION
# =============================================================================
DP = interp1d(np.asarray([0, 1225.0, 1675.0, 2125.0, 2575.0, 5000]), np.asarray([1.50823521,1.50823521,1.27057052,1.29932794,2.43052786,2.43052786]), kind = 'linear')     
# =============================================================================

# LOAD MODEL GRID
# =============================================================================
Xgrid =  np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
print('Model coordinates loaded.')
# =============================================================================

# Set-up time period and timesteps for downscaling:
# =============================================================================
yr = int(sys.argv[1]) # get year from the job array
print(type(yr))
years = np.array([yr])

# APPLY TEMPERATURE BIAS CORRECTION
# =============================================================================
print('Beginning Temperature Bias Correction')
for year in years:
    print('Bias correcting:', str(year),'Temperature')
    
    # Load downscaled temperature file:
    inT = Dataset(os.path.join(downscaled_inputs,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    T_array = inT.variables['Temperature'][:]
    sys.stdout.flush()

    # Set up file for bias-corrected outputs:
    Corrected_T = np.empty((T_array.shape))
    Corrected_T_file = os.path.join(output_path,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc')
    
    # Get DOY index for each timestep:
    DOY = np.linspace(1,(len(T_array)/8),len(T_array))
    
    # Loop through each timestep:
    for dt in range(0,len(T_array)):
        # Corrected T   = Downscaled T + Bias-correction factor (factor of time)
        Corrected_T[dt] = T_array[dt] + DT(DOY[dt])
    
    # Save downscaled and bias-corrected temperature as netcdf file:
    save_to_netcdf(Corrected_T, 'Temperature', Corrected_T_file, year, Xgrid, Ygrid,3)   

# APPLY ACCUMULATION BIAS CORRECTION
# =============================================================================
    print('Bias correcting:', str(year),'Precip')
    
    # Load downscaled precipitation file:
    inP = Dataset(os.path.join(downscaled_inputs,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    P_array = inP.variables['Precipitation'][:]
    sys.stdout.flush()
    
    # Load DEM:
    Zgrid =  np.loadtxt(os.path.join(Elev_inputs,'DEM_' + str(Glacier_ID) + '_' + str(year) + '.txt'))
    
    # Calculate correction factors based on elevation:
    Cz = DP(Zgrid)

    # Set up file for bias-corrected outputs:
    Corrected_P = np.empty((P_array.shape))
    Corrected_P_file = os.path.join(output_path,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc')

    # Partition rain vs snow based on downscaled temperature:    
    rainlocs = np.where(Corrected_T > R2S)
    
    # Correct all precipitation
    Corrected_P = P_array[:] * Cz
    
    # Replace rainy pixels with uncorrected precipitation
    Corrected_P[rainlocs] = P_array[rainlocs]
    
    # Save downscaled and bias-corrected precipitation as netcdf file:
    save_to_netcdf(Corrected_P, 'Precipitation', Corrected_P_file, year, Xgrid, Ygrid,3)   
   
    inT.close()
    inP.close()
    
print('Bias correction complete!! :)')

