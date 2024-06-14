# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 10:06:05 2023

Calculate 2007--2018 mass balance from the tuning outputs
DEMs used to calculate the 2007--2018 mass balance (Young et al. 2021) were from
(primarily) Sept 3 2007 and Oct 1 2018. 
Adjust variables start_year, end_year, DEMdate if needed.

We calculate the mass balance for this period and divide by the number of years (~11.09)
for each simulation.

The 2007-2018 modelled mass balance is saved in a txt file.

@author: katierobinson
"""

import numpy as np
from netCDF4 import Dataset
import os
import sys
import pandas as pd

# define the period over which the annual net mass balance will be calculated
start_year = 2007 
end_year = 2018

# path to model outputs
INPUT_PATH = '/home/krobin/scratch/ModelRuns/KRH/test_2024_06_14'

# glacier identifier (should be same as model outputs)
Glacier_ID = 'KRH'

# array that defines what parts of the grid to calculate mass balance for (off-glacier = NaN, glacier = finite)
glacier_grid = np.loadtxt('/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/InputGeometries/KRH_Tributaries.txt')

# Get sim number from job array
sim = int(sys.argv[1])
print('this is sim',sim)
print(type(sim))
sys.stdout.flush()

running_mb = np.zeros(glacier_grid.shape)
for year in range(start_year,end_year+1):
    print(year)
    sys.stdout.flush()
    dataset = Dataset(os.path.join(INPUT_PATH,'Netbalance_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
    mb = dataset.variables['Netbalance'][:]
    sys.stdout.flush()
    dataset.close()
    
    if year == start_year:
        dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
        DEMdate = np.where(dates == pd.Timestamp(str(year)+'-09-03T00'))[0][0]
        endofyear_mb = np.sum(mb[DEMdate:],axis=0)
        
    elif year == end_year:
        dates = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
        DEMdate = np.where(dates == pd.Timestamp(str(year)+'-10-01T00'))[0][0]
        endofyear_mb = np.sum(mb[:DEMdate],axis=0)
        
    else:
        endofyear_mb = np.sum(mb,axis=0)
        
    running_mb += endofyear_mb
    
num_years = len(pd.date_range(start= str(start_year) + '-09-03 00:00:00',end= str(end_year) + '-10-01 00:00:00',freq='D'))/365
annual_mb = running_mb/num_years
        
KW_mb = np.mean(annual_mb[np.isfinite(glacier_grid)])
print(KW_mb)
sys.stdout.flush() 

np.savetxt(os.path.join(INPUT_PATH,'sim' + str(sim) + '_mb.txt'),[KW_mb])




# Calculate modelled mass balance in each section of the Kaskawulsh Glacier bounded by the flux gates:

# Load map of flux gate sections:
#KRH_fluxgates = np.loadtxt('/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/InputGeometries/KRH_Fluxgates.txt')

#Zones = ['NA','CA','SW','SA','Main\nTrunk','KW0','KW1','KW2','KW3','KW4','KW5','KW\n1-2','KW\n4-5','All\nTs','Glacier-\nwide']
#bmod_KR = [np.round(np.mean(annual_mb[np.where(KRH_tributaries==4)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_tributaries==3)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_tributaries==2)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_tributaries==1)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_tributaries==5)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==0)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==9)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==8)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==7)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==6)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates==5)]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_fluxgates>=8)]),2),\
#        np.round(np.mean(annual_mb[np.where((KRH_fluxgates==6) | (KRH_fluxgates==5))]),2),\
#        np.round(np.mean(annual_mb[np.where(KRH_tributaries<=4)]),2),\
#        np.round(np.mean(annual_mb[np.isfinite(KRH_tributaries)]),2)]

#data = {'Zone': Zones, 'Bmod': bmod_KR}
#df = pd.DataFrame(data)
#df.to_csv(os.path.join(INPUT_PATH,'sim' + str(sim) + '_Bmod.txt'), sep='\t', index=False)
