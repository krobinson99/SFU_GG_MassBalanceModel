# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:52:00 2023

Script for downscaling NARR Temperature and Precipitation.
Based on downscaling routine from Young et al. (2021).

Set configuration options in downscaling_config.py

Inputs:
    NARR 3-hourly surface air temperature
    NARR daily surface precipitation
    DEM for NARR grid
    Coordinates for model grid (X,Y,Z)
Outputs:
    Downscaled T/P as netcdf files

@author: katierobinson
"""
      
# Import packages:
import numpy as np
from scipy import interpolate
import netCDF4
from pyproj import Proj
from netCDF4 import Dataset
import sys, os

# Import parameters from config file
from downscaling_config import Glacier_ID, UTM, NARR_subregions, time_step
from downscaling_config import Climate_inputs, Coarse_DEM_input, Easting_grid, Northing_grid, Elev_inputs, OUTPUT_PATH, Model_functions

# Import functions for model:
sys.path.insert(1,Model_functions)
from model_functions import write_config_file, save_to_netcdf, closest_node
from model_functions import precip_downscaling, T_downscale_funkfest

# Save configuration file for this run to output directory:
write_config_file(OUTPUT_PATH,"downscaling_config.py")

# Load (time-invariant) coarse NARR grid:
# =============================================================================
NARR_DEM = Dataset(Coarse_DEM_input, "r")
coarse_elev =  NARR_DEM.variables['hgt'][0,:,:] # surface geopotential height (time invariant), units = m a.s.l.

lons = NARR_DEM.variables['lon'][:]
lats = NARR_DEM.variables['lat'][:]
units = NARR_DEM.variables['time'].units
sys.stdout.flush()

Projection = Proj('+proj=utm +zone=' +str(UTM) + ' +ellps=WGS84', preserve_units=False)
UTMx, UTMy = Projection(lons, lats)  # converts lat/lons of coarse NARR grid to easting, northing on WGS84 projection.

#create list of array positions ungridded
UTMx_list = UTMx.ravel()
UTMy_list = UTMy.ravel()
print("Coarse NARR grid loaded")
# =============================================================================


# Load grid to downscale onto:
# =============================================================================
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
print('Model coordinates loaded.')
# =============================================================================


# Set-up time period and timesteps for downscaling:
# =============================================================================
yr = int(sys.argv[1]) # get year from the job array
print(type(yr))
years = np.array([yr])

time_steps = np.arange(0,24,time_step)
# =============================================================================

# Begin looping through years: 
# =============================================================================
for year in years:
    print('Starting downscaling for:',year)

    
    # Load inputs for downscaling (DEM, Temp, Precip, Geopotential Height) 
    # =========================================================================
    Zgrid = np.loadtxt(os.path.join(Elev_inputs,'DEM_' + str(Glacier_ID) + '_' + str(year) + '.txt'))
    nanlocs = np.where(np.isnan(Zgrid))
    print("Model DEM loaded.") 
    
    # Load NARR Geopotential Heights, Temperature, and Precip for year:
    File_elev_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_hgt.' + str(year) + '.nc')
    File_temp_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_air.' + str(year) + '.nc')
    File_precip_in = os.path.join(Climate_inputs,str(Glacier_ID) + '_apcp.' + str(year) + '.nc')
    
    # Geopotential Height:
    inH = Dataset(os.path.join(Climate_inputs,str(Glacier_ID) + '_hgt.' + str(year) + '.nc'), "r")
    sys.stdout.flush()
    H_array = inH.variables['hgt'][:]
    
    # Temperature:
    inT = Dataset(File_temp_in, "r")
    sys.stdout.flush()
    T_array = inT.variables['air'][:]
    
    # Precipitation:
    inP = Dataset(File_precip_in, "r")
    sys.stdout.flush()
    P_array = inP.variables['apcp'][:]
    print("NARR Geopotential, Temperature, Precipitation inputs loaded.") 
    sys.stdout.flush()
    # =========================================================================               
    
    
    # Prep output variables
    # ========================================================================= 
    Downscaled_T = np.empty((T_array.shape[0],Xgrid.shape[0],Xgrid.shape[1]))
    Downscaled_T_file = os.path.join(OUTPUT_PATH,'Temperature_' + str(Glacier_ID) + '_' + str(year) + '.nc')
    
    Downscaled_P = np.empty((T_array.shape[0],Xgrid.shape[0],Xgrid.shape[1]))
    Downscaled_P_file = os.path.join(OUTPUT_PATH,'Precipitation_' + str(Glacier_ID) + '_' + str(year) + '.nc')

    # Get timestamps for T (3 hourly) and P (daily)
    dt_daily = netCDF4.num2date(inP.variables['time'][:], units) # 365 timesteps (daily for 1 year)
    dt_3hourly = netCDF4.num2date(inH.variables['time'][:], units)   # 2920 timesteps (3 hourly for 1 year)
    date_list = []
    for timestep in dt_3hourly:
        date_list.append(timestep.replace(hour=0, minute=0, second=0, microsecond=0)) # sets hour to zero for all dates (still 2920 timesteps)
    # ========================================================================= 
    
    
    # Start looping through timesteps:
    # ========================================================================= 
    for date in dt_daily:
        print(date)
        sys.stdout.flush()
       
        # Get values for current iteration         
        hourly_indices = np.where(np.array(date_list) == date)
        daily_index = np.where(dt_daily == date)

        dailyH = H_array[hourly_indices] #Selects all pressure levels and gridcells for current day. Shape=(8,29,6,6)
        dailyT = T_array[hourly_indices] # (8,29,6,6)
        dailyP = P_array[daily_index][0]/1000  # (6,6)  # Convert daily precipitation from mm w.e. to m w.e.
    
# =============================================================================
#         PRECIP DOWNSCALING: 
# =============================================================================
        # Get coefficients from NARR Precip
        r_beta2, b_coeffs, b0 = precip_downscaling(coarse_elev.ravel()[NARR_subregions], dailyP.ravel()[NARR_subregions], UTMx_list[NARR_subregions], UTMy_list[NARR_subregions]) 
        # Calculate P_local across model grid using coeffs
        Plocal = (b0 + (b_coeffs[0] * Xgrid) + (b_coeffs[1] * Ygrid) + (b_coeffs[2] * (Xgrid * Ygrid)) + (b_coeffs[3] * (Xgrid**2)) + (b_coeffs[4] * (Ygrid**2)) + (b_coeffs[5] * Zgrid))*r_beta2 
        # Correct for negative precip values in areas where statistical model predicts less than zero value on mx + b regression curve
        Plocal[np.where(Plocal<0)] = 0       

        Pregional = np.empty(Plocal.shape)
        
         #Looping through each timestep
        for i in range(0, len(time_steps)):
            
            # Get Geopotential height and temperature vals for current timestep
            hourlyH = dailyH[i]      
            hourlyT = dailyT[i]
            
            # Apply T_downscale_funkfest function to get inversion tracking downscaled T functions for every reanalysis grid point
            xi_list, yi_list, xi_list_inver, yi_list_inver, funclist, funclist_inver, inversion_list, y0func, y0func_inver, Lfunc, Lfunc_inver = T_downscale_funkfest(hourlyT, hourlyH, UTMx_list, UTMy_list) 
            
            # Loop over every gridcell in downscaled reanalysis, get downscaled P and T
            for cell in range(0, len(np.where(np.isfinite(Zgrid))[0])):
                x = np.where(np.isfinite(Zgrid))[0][cell]
                y =  np.where(np.isfinite(Zgrid))[1][cell]
                
                #Get closest NARR grid point for appropriate downscaling T values
                downscaled_cell = np.asarray(([Xgrid[x,y]], [Ygrid[x,y]]))
                NARR_cell = closest_node(downscaled_cell, UTMx_list, UTMy_list) #Finds which NARR gridcell (36 total) is closest to the gridcell being downscaled.
            
                #use index to get nearest grid point in u, w notation
                u = int((np.where(UTMx == UTMx_list[NARR_cell])[0])[0])
                w = int((np.where(UTMy == UTMy_list[NARR_cell])[1])[0])

                
# =============================================================================
#                 TEMPERATURE DOWNSCALING
# =============================================================================
                if inversion_list[u][w] == 0:
                    # Interpolated lapse rate*elev + interpolated sea level temp
                    k = interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], Lfunc) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], y0func)                         
                        
                else:
                    
                    if Zgrid[x,y] < yi_list[u][w][0]:
                        k = interpolate.bisplev(Ygrid[x,y],Xgrid[x,y], Lfunc_inver) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], y0func_inver) 
                    else:
                        k = interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], Lfunc) * Zgrid[x,y] + interpolate.bisplev(Ygrid[x,y], Xgrid[x,y], y0func)  
                        
                K = (k - 273.15)
                Downscaled_T[hourly_indices[0][i],x,y] = K
                
                # Now get P_regional (Precip for nearest NARR node):
                if i == 0: #calculate for first timestep only
                    Pregional_cell = dailyP[u][w]*(1-r_beta2)
                    Pregional[x,y] = Pregional_cell
                else:
                    pass
            
            Downscaled_T[hourly_indices[0][i]][nanlocs] = np.nan
            
            if i == 0:
                Pdownscaled = (Plocal + Pregional)/len(time_steps) # Split daily precip throughout day
                Pdownscaled[nanlocs] = np.nan
            else:
                pass
            
            Downscaled_P[hourly_indices[0][i],:,:] = Pdownscaled
        
    # Save downscaled T/P as netcdf
    save_to_netcdf(Downscaled_P, 'Precipitation', Downscaled_P_file, year, Xgrid, Ygrid,3) 
    save_to_netcdf(Downscaled_T, 'Temperature', Downscaled_T_file, year, Xgrid, Ygrid,3)        
 
    inH.close()
    inP.close()
    inT.close()
 
print(Glacier_ID + 'Downscaling Complete')
    
