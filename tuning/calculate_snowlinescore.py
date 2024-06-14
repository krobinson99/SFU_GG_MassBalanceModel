# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 14:12:37 2023

Calculate snowline scores for each simulation by comparing modelled snow cover
with observed snow cover rasters delineated from satellite images of the Kaskawulsh Glacier.

Output: Simulation scores for each delineated satellite image:
- 1 score for the whole glacier
- 1 score for each tributary (SA, SW, NA, CA, TR)
- weight for each image (based on how much of the glacier area was visible in the satellite image)
- output file simN_snowlinescores.txt has shape (8, 53) = (8 scores, 53 images)

Final score is the time-weighted average of individual image scores (calculated in tuning_stageII.py)

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset
import sys
import os

# =============================================================================
# Config
# =============================================================================

path_to_obs_snowlines = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/Snowlines'   # directory with rasterized observed snowlines (.npy files)
MODEL_OUTPUTS = '/home/krobin/scratch/ModelRuns/KRH/test_2024_06_14'         # directory of model outputs (this is also where S.S. will be saved)
Glacier_ID = 'KRH'                                                                                          # glacier identifier (should be same as model outputs)

# array that defines the tributaries of the Kaskawulsh Glacier (integer from 1-5 for each tributary)
KRH_tributaries = np.loadtxt('/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/InputGeometries/KRH_Tributaries.txt')
nanlocs = np.where(np.isnan(KRH_tributaries))

snowline_years = np.arange(2012,2019+1)  # satellite images from 2012-2019 were delineated

# Get sim number from job array
sim_id = int(sys.argv[1])
print('this is sim',sim_id)
print(type(sim_id))
sims = np.arange(sim_id,sim_id+1) 

# =============================================================================
# Functions
# =============================================================================

# Get array of dates where there is an observed snowline:
all_snow_depth_dates = []
for file in os.listdir(path_to_obs_snowlines):
    if file.endswith('.npy'):
        date = pd.Timestamp(year=int(file[15:19]), month=int(file[20:22]), day=int(file[23:25]))
        all_snow_depth_dates.append(date)

def plotmodelsnowline(timestamp):
    snow = snow_depth[np.where(snow_depth_dates==timestamp)[0][0]]
    plt.figure(figsize=(7,5))
    plt.contourf(Xgrid,np.flipud(Ygrid),snow,cmap='BuPu',levels=np.linspace(0,1.5,16))
    plt.axis('equal')
    legend = plt.colorbar()
    legend.ax.set_ylabel('Accumulation (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    legend.ax.tick_params(labelsize=14)
    plt.title('Modelled snow depth: ' + str(timestamp),fontsize=14)
    plt.tight_layout()


def calculate_snowline_score(obs_snowline,model_snowdepth):
    '''
    calculate score and weight for each delineated satellite image
    
    score: percentage of matching cells out of all the cells where sfc condition is
    known (ice or snow, does not count transition zone)
    
    weight: the fraction of cells out of the total grid that have info
    '''
    # Set cells with no info == 1
    #snow = 2
    #ice = 0
    #transition zone = 1
    #no info = nan
    obs_snowline[np.where(obs_snowline == -0.5)] = np.nan
    
    # Convert model snowdepth to a snowline
    model_snowdepth[np.where(model_snowdepth > 0)] = 2
    
    # get the absolute difference between modelled and observed snowline:
    # 0 = match
    # 1 = transition zone
    # 2 = mismatch (ice where snow should be, or vice versa)
    difference = np.abs(obs_snowline - model_snowdepth)
    norm = np.sqrt(np.nansum(np.square(difference)))
    
    # calculate score: percentage of matching cells (==0)
    num_matching_cells = len(np.where(difference==0)[0])
    total_cells_snowandice = len(np.where(obs_snowline==0)[0]) + len(np.where(obs_snowline==2)[0])
    image_score = num_matching_cells/total_cells_snowandice
    
    # calculate weight: the fraction of cells out of the total grid that have info
    weight = total_cells_snowandice/np.where(np.isfinite(KRH_tributaries))[0].shape[0]

    
    # calculate individual trib scores: (percent of each trib with info)
    tribscores = []
    for trib in range(1,6):
        
        matchingcells_intrib = (len(np.where(difference[np.where(KRH_tributaries==trib)] == 0)[0]))
        
        cellsintrib_with_info = (len(np.where(obs_snowline[np.where(KRH_tributaries==trib)] == 0)[0]) + len(np.where(obs_snowline[np.where(KRH_tributaries==trib)] == 2)[0]))
        
        total_cellsintrib = len(np.where(KRH_tributaries==trib)[0])
        
        if cellsintrib_with_info/total_cellsintrib <= 0.2:
            tribscore = np.nan
        else:
            tribscore = matchingcells_intrib/cellsintrib_with_info
        tribscores.append(tribscore)
        
    SAscore, SWscore, CAscore, NAscore, TRscore = tribscores
    
    return image_score, weight, SAscore, SWscore, CAscore, NAscore, TRscore, norm
    
# =============================================================================
# Calculate snowline scores
# ============================================================================= 
    
# Get snow depth array for each year in each simulation:
for sim in sims:
    print('Simulation: ',sim)
    simscores = np.zeros((8,len(all_snow_depth_dates))) # (image score (whole glacier), weight, South Arm score, )
    simscores[:] = np.nan
    
    for year in snowline_years[:-1]:
        print('Accumulation season starting in September',year)
        snow_depth_dates = pd.date_range(start=str(year)+'-09-01 00:00:00',end=str(year+1)+'-08-31 21:00:00',freq='D')
        
        dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
        dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq='D')
        
        inMB1 = Dataset(os.path.join(MODEL_OUTPUTS,'Snowdepth' + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
        inMB2 = Dataset(os.path.join(MODEL_OUTPUTS,'Snowdepth' + '_' + str(Glacier_ID) + '_' + str(year+1) + '_' + str(sim) + '.nc'),'r')
        
        mb_array_year1 = inMB1.variables['Snowdepth'][:]
        sys.stdout.flush()
        inMB1.close()
    
        mb_array_year2 = inMB2.variables['Snowdepth'][:]
        sys.stdout.flush()
        inMB2.close()
        
        mb_yr1 = mb_array_year1[np.where(dates_yr1 == pd.Timestamp(str(year)+'-09-01T00'))[0][0]:]
        mb_yr2 = mb_array_year2[:np.where(dates_yr2 == pd.Timestamp(str(year+1)+'-09-01T00'))[0][0]]
        mb_hydroyear = np.concatenate((mb_yr1,mb_yr2),axis=0)

        snowdepth_change = np.zeros(mb_hydroyear.shape)
        snow_depth = np.zeros(mb_hydroyear.shape)
        for i in range(0,len(mb_hydroyear)-1):
            snowdepth_change[i+1] += mb_hydroyear[i+1] - mb_hydroyear[i]
            snow_depth[i+1] = snow_depth[i] + snowdepth_change[i+1]
            snow_depth[i+1][np.where(snow_depth[i+1] < 0)] = 0
            
        #Now loop through all the snowlines in this year and get the individual image score 
        # Get snowlines between sept 1 of current year and aug 31 of following year
        for file in os.listdir(path_to_obs_snowlines):
            if file.endswith('.npy'):
                date = pd.Timestamp(year=int(file[15:19]), month=int(file[20:22]), day=int(file[23:25]))
              
                acc_season_start = pd.Timestamp(year=year, month=9, day=1)
                acc_season_end = pd.Timestamp(year=year+1, month=8, day=31)
                
                if date >= acc_season_start:
                    if date <= acc_season_end:
                        print(file,date)
                        
                        obs_snowline = np.load(os.path.join(path_to_obs_snowlines,file))
                        model_snowdepth = snow_depth[np.where(snow_depth_dates==date)[0][0]]
                        model_snowdepth[nanlocs] = np.nan
                        
                        image_score, weight, SAscore, SWscore, CAscore, NAscore, TRscore, norm = calculate_snowline_score(obs_snowline,model_snowdepth)
                        dateindex = np.where(np.array(all_snow_depth_dates)==date)[0][0]
                        
                        simscores[0][dateindex] = image_score
                        simscores[1][dateindex] = weight
                        simscores[2][dateindex] = SAscore
                        simscores[3][dateindex] = SWscore
                        simscores[4][dateindex] = CAscore
                        simscores[5][dateindex] = NAscore
                        simscores[6][dateindex] = TRscore
                        simscores[7][dateindex] = norm
                    
                    else:
                        pass
                else:
                    pass
          
    np.savetxt(os.path.join(MODEL_OUTPUTS,'sim' + str(sim) + '_snowlinescores.txt'),simscores)        
