# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:02:43 2021

Distributed glacier mass-balance model originally 
developped by E. Young (Young et al. 2021), updated 
by K. Robinson (Robinson et al. 2024).

Set configuration options in mbm_config.py

Inputs:
    Melt model parameters (MF, asnow, aice)
    Downscaled/bias-corrected NARR 3-hourly surface air temperature
    Downscaled/bias-corrected NARR 3-hourly surface precipitation
    3 hourly potential direct incoming solar radiation
    Coordinates for model grid (X,Y)
    Surface Type Grid (1 = off-glacier, 0 = on-glacier, NaN = not in the domain.)
    Debris cover grid

Outputs:
    Distributed runoff, refreezing, and mass balance as NetCDF files,
    (1 file per variable per year with daily timesteps)

@author: katierobinson
"""

# Import packages:
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import sys
import os

# Import parameters from config file
from mbm_config import start_year, end_year, timestep, params_file, R2S, Glacier_ID, \
Model_functions, Precip_inputs, Temp_inputs, Solar_inputs, Easting_grid, Northing_grid, \
Sfc_grid, debris_parameterization, debris_map, cleanice_melt, peak_melt, peak_melt_thickness, \
transition_thickness, b0, k, OUTPUT_PATH, output_runoff, output_refreezing, output_accumulation, output_massbalance

# Import functions for model:
sys.path.insert(1,Model_functions)
from model_functions import write_config_file, save_to_netcdf
from model_functions import debris, Calculate_Pmean, max_superimposed_ice,  \
max_superimposed_ice_finalyear, rain_refreezing, updated_superimposed_ice, MassBalance

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"mbm_config.py")

# Get sim param from the job array
# =============================================================================
sim = int(sys.argv[1])
print('this is sim',sim)
print(type(sim))

# Load model grid:
# =============================================================================
years = np.arange(start_year,end_year+1)
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Sfc = np.loadtxt(Sfc_grid)
print('Model coordinates loaded.')
sys.stdout.flush()
# =============================================================================

# Load melt parameters for this run:
# =============================================================================
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice, asnow, MF = params[sim][1:]
# Save params for this run in a text file just in case:
np.savetxt(os.path.join(OUTPUT_PATH,'sim' + str(sim) + '_params.txt'),params[sim][1:])
print('Parameters for this run:\naice = ',aice,'\nasnow = ',asnow,'\nMF = ',MF)
sys.stdout.flush()
# =============================================================================
  
# Set up the debris parameterization:
# =============================================================================
debris_m = debris(debris_parameterization,debris_map,Sfc,cleanice_melt,peak_melt,peak_melt_thickness,transition_thickness,b0,k)
print('Debris parameterization loaded with option:',debris_parameterization)
sys.stdout.flush()
# =============================================================================

# Calculate mean annual precipitation (for refreezing):
# =============================================================================
Pmean = Calculate_Pmean(years,Glacier_ID,Precip_inputs,Temp_inputs,Sfc)
print('Mean annual snowpack calculated') 
sys.stdout.flush()
# =============================================================================

for year in years:
    print('Running model for',year)
    sys.stdout.flush()
    dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq=str(timestep)+'H')
    
    # Load inputs for model (Temperature, Precip, Solar Radiation) 
    # =========================================================================
    
    # Temperature:
    inT = Dataset(os.path.join(Temp_inputs,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    T_array = inT.variables['Temperature'][:]
    sys.stdout.flush()
    
    # Precipitation:
    inP = Dataset(os.path.join(Precip_inputs,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    P_array = inP.variables['Precipitation'][:]
    sys.stdout.flush()
    
    # Solar radiation:
    inS = Dataset(os.path.join(Solar_inputs,'Solar_' + str(Glacier_ID) + '_' + str(year) + '.nc'),'r')
    S_array = inS.variables['SolarRadiation'][:]
    sys.stdout.flush()
    
    print(year,'inputs (temperature, precipitation, solar) loaded.') 
    sys.stdout.flush()
    # =========================================================================     

    # Get maximum amount of superimposed ice that can form each hydrologic year:
    # =========================================================================
    if year == years[-1]:
        SImax = max_superimposed_ice_finalyear(year, P_array, T_array, timestep, Pmean)
    else:
        SImax = max_superimposed_ice(year, P_array, T_array, timestep, Glacier_ID, Pmean, Precip_inputs, Temp_inputs)

    # Set up output files and trackers for snowpack, potential superimposed ice:
    # =========================================================================
    IceMelt = np.zeros(T_array.shape)
    SnowMelt = np.zeros(T_array.shape)
    RefrozenMelt = np.zeros(T_array.shape)
    RefrozenRain = np.zeros(T_array.shape)
    MassBal = np.zeros(T_array.shape)
    
    if year == years[0]:
        # Trackers should have +1 extra timestep to account for carry over values for following year:
        Snowpack_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        PotentialSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        CurrentSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
    else:
        # Save final snowpack and leftover potential superimposed ice from previous year here.
        Snowpack_carryover = np.array(Snowpack_tracker[-1])
        PotentialSI_carryover = np.array(PotentialSI_tracker[-1])
        CurrentSI_carryover = np.array(CurrentSI_tracker[-1])
        
        # Reset snowpack for rest of the year to zero
        Snowpack_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        PotentialSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        CurrentSI_tracker = np.zeros((T_array.shape[0]+1,T_array.shape[1],T_array.shape[2]))
        
        # Add carry over values from previous year to beginning of tracker
        Snowpack_tracker[0] = Snowpack_carryover
        PotentialSI_tracker[0] = PotentialSI_carryover
        CurrentSI_tracker[0] = CurrentSI_carryover

    # Set up timestepping (loop through every timestep in year)
    # =========================================================================
    for timestamp in range(0,len(dates)):
        print(dates[timestamp])
        sys.stdout.flush()
        
        # Get current snowpack, potential SI, and actual SI volumes for this timestep:
        # =====================================================================
        # Renew cold content at beginning of hydrological year:
        if dates[timestamp] == pd.Timestamp(str(year)+'-10-01T00'):
            PotentialSI_tracker[timestamp] += SImax
        else:
            pass
        
        New_snowfall =  np.array(P_array[timestamp,:,:])
        New_snowfall[np.where(T_array[timestamp,:,:] > R2S)] = 0
        Snowpack_tracker[timestamp,:,:] += New_snowfall
        

        # Calculate Melt:
        # =====================================================================        
        Msnow, Mice, Refreezing, SI_out, SP_out = MassBalance(MF,asnow,aice,T_array[timestamp,:,:],S_array[timestamp,:,:],Snowpack_tracker[timestamp,:,:],PotentialSI_tracker[timestamp,:,:],debris_m,debris_parameterization,Sfc,CurrentSI_tracker[timestamp,:,:])
        
        # Calculate refreezing of rain (if any)
        # =====================================================================  
        refrozen_rain = rain_refreezing(P_array[timestamp,:,:],T_array[timestamp,:,:],R2S,SI_out,SP_out)
        SI_out -= refrozen_rain # Subtract the amount of rain that is refrozen from the potentialSI

        # Update output arrays for this timestep:
        # Total Melt = Snow Melt + Ice Melt
        # Net ablation = total melt - refreezing
        # Net balance = Accumulation (snow + rain that freezes) - net ablation (snowmelt + icemelt - refrozen snowmelt)
        # ===================================================================== 
        IceMelt[timestamp,:,:] = Mice
        SnowMelt[timestamp,:,:] = Msnow
        RefrozenMelt[timestamp,:,:] = Refreezing
        RefrozenRain[timestamp,:,:] = refrozen_rain
        MassBal[timestamp,:,:] = New_snowfall + RefrozenRain[timestamp,:,:] - ((Msnow - RefrozenMelt[timestamp,:,:]) + Mice)  

        # Update snowpack and superimposed ice trackers for next timestep:
        # ===================================================================== 
        Snowpack_tracker[timestamp+1,:,:] = SP_out
        PotentialSI_tracker[timestamp+1,:,:] = SI_out   
        CurrentSI_tracker[timestamp+1,:,:] = updated_superimposed_ice((RefrozenMelt[timestamp,:,:]+RefrozenRain[timestamp,:,:]),IceMelt[timestamp,:,:],CurrentSI_tracker[timestamp,:,:])
        
    # Save outputs for the year before starting next year:
    # =========================================================================
    print('Saving model outputs for',year)
    sys.stdout.flush()

    # Get daily values for each output:
    
    # Total snowmelt:
    SnowMelt_daily = np.array([sum(SnowMelt[i:i+8]) for i in range(0, len(SnowMelt), 8)])
    
    # Refreezing of snowmelt:
    RefrozenMelt_daily = np.array([sum(RefrozenMelt[i:i+8]) for i in range(0, len(RefrozenMelt), 8)])
    
    # Net snowmelt (snow melt that runs off = total snowmelt - refreezing)
    NetSnowMelt = np.subtract(SnowMelt,RefrozenMelt)
    NetSnowMelt_daily = np.array([sum(NetSnowMelt[i:i+8]) for i in range(0, len(NetSnowMelt), 8)])
    
    # Glacier ice melt:
    Glacier_IceMelt = np.subtract(IceMelt,CurrentSI_tracker[:-1])
    Glacier_IceMelt[np.where(Glacier_IceMelt < 0)] = 0
    Glacier_IceMelt_daily = np.array([sum(Glacier_IceMelt[i:i+8]) for i in range(0, len(Glacier_IceMelt), 8)])
    
    # Superimposed ice melt:
    Superimposed_IceMelt = np.subtract(IceMelt,Glacier_IceMelt)
    Superimposed_IceMelt_daily = np.array([sum(Superimposed_IceMelt[i:i+8]) for i in range(0, len(Superimposed_IceMelt), 8)])
            
    # Rainfall:
    Rain = np.array(P_array)
    Rain[np.where(T_array <= R2S)] = 0
    Rain_daily = np.array([sum(Rain[i:i+8]) for i in range(0, len(Rain), 8)])
    
    # Refreezing of rain:
    RefrozenRain_daily = np.array([sum(RefrozenRain[i:i+8]) for i in range(0, len(RefrozenRain), 8)])
    
    # Net rainfall (rain that runs off = total rainfall - refreezing)
    RainRunoff = np.subtract(Rain,RefrozenRain)
    RainRunoff_daily = np.array([sum(RainRunoff[i:i+8]) for i in range(0, len(RainRunoff), 8)])
    
    # Accumulation:
    Accumulation = np.array(P_array)
    Accumulation[np.where(T_array > R2S)] = 0
    Accumulation_daily = np.array([sum(Accumulation[i:i+8]) for i in range(0, len(Accumulation), 8)])
    
    # Mass balance:
    MassBal_daily = np.array([sum(MassBal[i:i+8]) for i in range(0, len(MassBal), 8)])

    # For the tracker arrays, just need to get the value once a day:
    Snowpack = np.array([sum(Snowpack_tracker[:-1][i:i+1]) for i in range(0, len(Snowpack_tracker[:-1]), 8)])
    Superimposedice = np.array([sum(CurrentSI_tracker[:-1][i:i+1]) for i in range(0, len(CurrentSI_tracker[:-1]), 8)])
    Ptau = np.array([sum(PotentialSI_tracker[:-1][i:i+1]) for i in range(0, len(PotentialSI_tracker[:-1]), 8)])
    
    # Save to netcdf: 
    if output_runoff == True:
        save_to_netcdf(NetSnowMelt_daily, 'Netsnowmelt', os.path.join(OUTPUT_PATH,'Netsnowmelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(Glacier_IceMelt_daily, 'Glaciericemelt', os.path.join(OUTPUT_PATH,'Glaciericemelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(RainRunoff_daily, 'Rainrunoff', os.path.join(OUTPUT_PATH,'Rainrunoff_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(Superimposed_IceMelt_daily, 'Superimposedicemelt', os.path.join(OUTPUT_PATH,'Superimposedicemelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
    else:
        pass

    if output_refreezing == True:
        save_to_netcdf(SnowMelt_daily, 'Snowmelt', os.path.join(OUTPUT_PATH,'Snowmelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(RefrozenMelt_daily, 'Refrozenmelt', os.path.join(OUTPUT_PATH,'Refrozenmelt_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(Rain_daily, 'Rain', os.path.join(OUTPUT_PATH,'Rain_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(RefrozenRain_daily, 'Refrozenrain', os.path.join(OUTPUT_PATH,'Refrozenrain_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24) 
        save_to_netcdf(Ptau, 'Ptau', os.path.join(OUTPUT_PATH,'Ptau_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(Superimposedice, 'Superimposedice', os.path.join(OUTPUT_PATH,'Superimposedice_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)       
    else:
        pass

    if output_accumulation == True:
        save_to_netcdf(Accumulation_daily, 'Accumulation', os.path.join(OUTPUT_PATH,'Accumulation_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
        save_to_netcdf(Snowpack, 'Snowdepth', os.path.join(OUTPUT_PATH,'Snowdepth_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
    else:
        pass

    if output_massbalance == True:
        save_to_netcdf(MassBal_daily, 'Netbalance', os.path.join(OUTPUT_PATH,'Netbalance_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'), year, Xgrid, Ygrid,24)
    else:
        pass

    inT.close()
    inP.close()
    inS.close()

print('RUN COMPLETE!')
