# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 10:07:22 2024

Plotting script for JoG Manuscript 1

@author: katierobinson
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from plotting_functions import load_daily_timeseries, load_distributed_fields, \
shiftedColorMap, distributed_mass_balance, hydrograph, distributed_runoff, massbalance_timeseries_average, \
cumulative_and_annual_massbal, distributed_var

# =============================================================================
# Set paths to files needed for plotting:
# =============================================================================

reference_model_outputs = '/Users/katierobinson/Desktop/MassBalanceModel/Demo'     # parent directory where formatted timeseries and distributed fields are saved
Glacier_ID = 'KRH'                                                 # identifier for the glacier. All inputs/outputs will also have the same ID.

# time period:
start_year = 1979   # timeseries starts Oct 1 of start_year
end_year = 2022     # timeseries ends Sep 30 of end_year

# Input geometry
KW_grid_file = '/Users/katierobinson/Desktop/MassBalanceModel/Data/input_geometry/KRH/KRH_Tributaries.txt' # text file where Kaskawulsh Glacier is finite, otherwise NaN 
Easting_grid = '/Users/katierobinson/Desktop/MassBalanceModel/Data/input_geometry/KRH/KRH_Xgrid.txt' # Paths to text files defining Easting/Northing coords of every model gridcell
Northing_grid = '/Users/katierobinson/Desktop/MassBalanceModel/Data/input_geometry/KRH/KRH_Ygrid.txt'
Elevation_grid = '/Users/katierobinson/Desktop/MassBalanceModel/Data/input_geometry/KRH/Zgrids/DEM_KRH_2022.txt'
Sfc_grid = '/Users/katierobinson/Desktop/MassBalanceModel/Data/input_geometry/KRH/KRH_SfcType.txt'      # Path to text file where 1 = off-glacier, 0 = on-glacier, NaN = not in the domain.

# =============================================================================
# Load coordinates and domain files for plotting
# =============================================================================

# Input geometry
Xgrid = np.loadtxt(Easting_grid)
Ygrid = np.loadtxt(Northing_grid)
Zgrid = np.loadtxt(Elevation_grid)
Sfc = np.loadtxt(Sfc_grid)
KW_grid = np.loadtxt(KW_grid_file)

# time period
years = np.arange(start_year,end_year)

# catchment outline for distributed plots
Catchmentoutline = np.array(Sfc)
Catchmentoutline[np.where(np.isfinite(Sfc))] = 0
Catchmentoutline[np.where(np.isnan(Sfc))] = 1

All_glacierized_area = np.array(Sfc)
All_glacierized_area[np.where(Sfc==1)] = np.nan

# mass balance model output variables to be loaded:
varnames = ['massbal','totalsnowmelt','refreezing','netsnowmelt','glaciericemelt','superimposedice_melt','rain','refrozen_rain','rain_runoff','accumulation','snowdepth','Ptau','SI']

# =============================================================================
# Load the timeseries and distributed outputs formatted as .txt files from format_mbm_outputs.py
# =============================================================================

# catchment-wide timeseries:
# outputs have shape (N years, 366 days)
massbal_krh_ref, totalsnowmelt_krh_ref, refreezing_krh_ref, netsnowmelt_krh_ref, gl_icemelt_krh_ref, superimp_icemelt_krh_ref, rain_krh_ref, refrozen_rain_krh_ref, rain_runoff_krh_ref, accumulation_krh_ref, \
snowdepth_krh_ref, Ptau_krh_ref, SI_krh_ref = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'krh',years)

# standard deviation for the catchment-wide timeseries:
massbal_krh_ref_std, totalsnowmelt_krh_ref_std, refreezing_krh_ref_std, netsnowmelt_krh_ref_std, gl_icemelt_krh_ref_std, superimp_icemelt_krh_ref_std, rain_krh_ref_std, refrozen_rain_krh_ref_std, rain_runoff_krh_ref_std, accumulation_krh_ref_std, \
snowdepth_krh_ref_std, Ptau_krh_ref_std, SI_krh_ref_std = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'krh_std',years)

# Kaskawulsh-wide timeseries
massbal_kw_ref, totalsnowmelt_kw_ref, refreezing_kw_ref, netsnowmelt_kw_ref, gl_icemelt_kw_ref, superimp_icemelt_kw_ref, rain_kw_ref, refrozen_rain_kw_ref, rain_runoff_kw_ref, accumulation_kw_ref, \
snowdepth_kw_ref, Ptau_kw_ref, SI_kw_ref = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'kw',years)

# standard deviation for the Kaskawulsh-wide timeseries:
massbal_kw_ref_std, totalsnowmelt_kw_ref_std, refreezing_kw_ref_std, netsnowmelt_kw_ref_std, gl_icemelt_kw_ref_std, superimp_icemelt_kw_ref_std, rain_kw_ref_std, refrozen_rain_kw_ref_std, rain_runoff_kw_ref_std, accumulation_kw_ref_std, \
snowdepth_kw_ref_std, Ptau_kw_ref_std, SI_kw_ref_std = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'kw_std',years)

# glacier-wide timeseries
massbal_allgl_ref, totalsnowmelt_allgl_ref, refreezing_allgl_ref, netsnowmelt_allgl_ref, gl_icemelt_allgl_ref, superimp_icemelt_allgl_ref, rain_allgl_ref, refrozen_rain_allgl_ref, rain_runoff_allgl_ref, accumulation_allgl_ref, \
snowdepth_allgl_ref, Ptau_allgl_ref, SI_allgl_ref = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'allgl',years)

# standard deviation for the glacier-wide timeseries:
massbal_allgl_ref_std, totalsnowmelt_allgl_ref_std, refreezing_allgl_ref_std, netsnowmelt_allgl_ref_std, gl_icemelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, rain_allgl_ref_std, refrozen_rain_allgl_ref_std, rain_runoff_allgl_ref_std, accumulation_allgl_ref_std, \
snowdepth_allgl_ref_std, Ptau_allgl_ref_std, SI_allgl_ref_std = load_daily_timeseries(reference_model_outputs,'REF_MODEL',varnames,'allgl_std',years)

# distributed fields:
# outputs have shape (43 years x 230 x 329)
massbal_dist_ref, totalsnowmelt_dist_ref, refreezing_dist_ref, netsnowmelt_dist_ref, gl_icemelt_dist_ref, superimp_icemelt_dist_ref, rain_dist_ref, refrozen_rain_dist_ref, rain_runoff_dist_ref, accumulation_dist_ref, \
snowdepth_dist_ref, Ptau_dist_ref, SI_dist_ref = load_distributed_fields(reference_model_outputs,'REF_MODEL',varnames,years)

# =============================================================================
# Plotting reference model figures (figs. 5.17, 5.21, 5.22 from KR's thesis)
# =============================================================================

# Define custom colormaps for figures
shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.8298755186721992,stop=1,name='massbal') # custom colormap for mass balance field

# Fig 5.17 from KR thesis (distributed mass balance, 1980-2022)
distributed_mass_balance(np.arange(1980,2022),years,massbal_dist_ref,np.linspace(-10,2.05,30),Xgrid,Ygrid,Sfc,Catchmentoutline,KW_grid,reference_model_outputs)
distributed_var(np.arange(1980,2022),years,accumulation_dist_ref,np.arange(0,2.1,0.2),'Accumulation (m w.e. a$^{-1}$)','BuPu',Xgrid,Ygrid,Sfc,Catchmentoutline)
distributed_var(np.arange(1980,2022),years,refreezing_dist_ref,np.arange(0,0.151,0.02),'Refreezing (m w.e. a$^{-1}$)','RdPu',Xgrid,Ygrid,Sfc,Catchmentoutline)

# fig 5.21 from KR thesis (catchement-wide average hydrograph, 1980-2022)
catchment_runoff = hydrograph('m3','Catchment-wide average runoff: 1980-2022 ',np.arange(1980,2022),years,450,2.7,gl_icemelt_krh_ref, netsnowmelt_krh_ref, rain_runoff_krh_ref, superimp_icemelt_krh_ref,gl_icemelt_krh_ref_std, netsnowmelt_krh_ref_std, rain_runoff_krh_ref_std, superimp_icemelt_krh_ref_std,Sfc,reference_model_outputs) 

#fig 5.21 modified to show runoff averaged over the Kaskawulsh Glacier only:
kaskawulsh_runoff = hydrograph('m3','Kaskawulsh Glacier average runoff: 1980-2022 ',np.arange(1980,2022),years,450,2.7,gl_icemelt_kw_ref, netsnowmelt_kw_ref, rain_runoff_kw_ref, superimp_icemelt_kw_ref,gl_icemelt_kw_ref_std, netsnowmelt_kw_ref_std, rain_runoff_kw_ref_std, superimp_icemelt_kw_ref_std,KW_grid,reference_model_outputs) 

# fig 5.22 from KR thesis (distributed runoff components, 1980-2022)
distributed_runoff(np.arange(1980,2022),years,Sfc,Xgrid,Ygrid,Catchmentoutline,gl_icemelt_dist_ref,netsnowmelt_dist_ref,rain_runoff_dist_ref,superimp_icemelt_dist_ref,reference_model_outputs)

# cumulative annual balance curve
massbalance_timeseries_average('Glacierized area mass balance: ',np.arange(1980,2022),years,0.035,1, accumulation_allgl_ref, refrozen_rain_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref, gl_icemelt_allgl_ref, accumulation_allgl_ref_std, refrozen_rain_allgl_ref_std, netsnowmelt_allgl_ref_std, superimp_icemelt_allgl_ref_std, gl_icemelt_allgl_ref_std, refreezing_allgl_ref, All_glacierized_area)
massbalance_timeseries_average('Kaskawulsh mass balance: ',np.arange(1980,2022),years,0.035,1, accumulation_kw_ref, refrozen_rain_kw_ref, netsnowmelt_kw_ref, superimp_icemelt_kw_ref, gl_icemelt_kw_ref, accumulation_kw_ref_std, refrozen_rain_kw_ref_std, netsnowmelt_kw_ref_std, superimp_icemelt_kw_ref_std, gl_icemelt_kw_ref_std,refreezing_kw_ref,KW_grid)

# cumulative mass loss
cumulative_and_annual_massbal('Glacierized Area Mass Balance\n',years,np.arange(1980,2022),0.07,22,2.5,accumulation_allgl_ref, refrozen_rain_allgl_ref, gl_icemelt_allgl_ref, netsnowmelt_allgl_ref, superimp_icemelt_allgl_ref,All_glacierized_area) #plt.savefig('D:\Model Runs\REF_MODEL\Plots\Refmodel_MB_timeseries_1979-2022.pdf',bbox_inches='tight')
