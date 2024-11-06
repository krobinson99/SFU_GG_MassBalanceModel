# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:27:51 2024

@author: katierobinson
"""

import numpy as np
import pandas as pd
import os
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import cmocean
import gc
import heapq
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec

# =============================================================================
# Functions for formatting the mbm outputs
# =============================================================================

def load_hydrologic_year(sim,year,varname,model_outputs,Glacier_ID,stddev):
    '''
    Returns a given variable from Oct 1 -- Sept 30 (hydrological year)
    '''
        
    # Concatenate variable from Oct 1 of current year to Sept 30 of following year    
    dates_yr1 = pd.date_range(start= str(year) + '-01-01 00:00:00',end= str(year) + '-12-31 21:00:00',freq='D')
    dates_yr2 = pd.date_range(start= str(year+1) + '-01-01 00:00:00',end= str(year+1) + '-12-31 21:00:00',freq='D')
    inMB1 = Dataset(os.path.join(model_outputs,varname + '_' + str(Glacier_ID) + '_' + str(year) + '_' + str(sim) + '.nc'),'r')
    inMB2 = Dataset(os.path.join(model_outputs,varname + '_' + str(Glacier_ID) + '_' + str(year+1) + '_' + str(sim) + '.nc'),'r')

    if stddev == True:
        varname = str(varname[:-4])
        
    mb_array_year1 = inMB1.variables[varname][:]
    sys.stdout.flush()
    inMB1.close()
    
    mb_array_year2 = inMB2.variables[varname][:]
    sys.stdout.flush()
    inMB2.close()
    
    mb_yr1 = mb_array_year1[np.where(dates_yr1 == pd.Timestamp(str(year)+'-10-01T00'))[0][0]:]
    mb_yr2 = mb_array_year2[:np.where(dates_yr2 == pd.Timestamp(str(year+1)+'-10-01T00'))[0][0]]
    mb_hydroyear = np.concatenate((mb_yr1,mb_yr2),axis=0)
    
    return mb_hydroyear


def process_model_timeseries(model_outputs, glacier_ID, sim, years, glacier_grid, std_dev):
    '''
    model_outputs = directory with the model outputs
    glacier_ID = identifier for the glacie (string)
    sim = integer corresponding to the simulation number that is being processed
    years = years to be processed (array)
    glacier_grid = array where gridcells to be averaged over are finite, other cells are NaN
    std_dev = boolean, True if calculating standard deviation, False if just getting regular variable
    
    returns (all units = m w.e. day^-1):
    daily runoff components:
        net snowmelt (total snow melt minus refreezing)
        superimposed ice melt
        glacier ice melt
        rain runoff (total rain minus refrozen rain)
    
    daily refreezing components:
        total snow melt
        refrozen snow melt
        total rainfall
        refrozen rainfall
        potential retention mass (Ptau)
        superimposed ice thickness
        
    daily mass balance:
        net mass balance
        accumulation
        
    other:
        snowpack thickness
    '''
    
    massbal_l = []
    snowmelt_l = []
    refrozenmelt_l = []
    netsnowmelt_l = []
    glaciermelt_l = []
    SImelt_l = []
    rain_l = []
    refrozenrain_l = []
    rainrunoff_l = []
    accumulation_l = []
    snowdepth_l = []
    Ptau_l = []
    SI_l = []
    
    if std_dev == True:
        var_suffix = '_std'
    else:
        var_suffix = ''
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim, year, 'Netbalance' + var_suffix, model_outputs, glacier_ID, std_dev)
        snowmelt = load_hydrologic_year(sim, year, 'Snowmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        refrozenmelt = load_hydrologic_year(sim, year, 'Refrozenmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        netsnowmelt = load_hydrologic_year(sim, year, 'Netsnowmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        glaciermelt = load_hydrologic_year(sim, year, 'Glaciericemelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        SImelt = load_hydrologic_year(sim, year, 'Superimposedicemelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        rain = load_hydrologic_year(sim, year, 'Rain' + var_suffix, model_outputs, glacier_ID, std_dev)
        refrozenrain = load_hydrologic_year(sim, year, 'Refrozenrain' + var_suffix, model_outputs, glacier_ID, std_dev)
        rainrunoff = load_hydrologic_year(sim, year, 'Rainrunoff' + var_suffix, model_outputs, glacier_ID, std_dev)
        accumulation = load_hydrologic_year(sim, year, 'Accumulation' + var_suffix, model_outputs, glacier_ID, std_dev)
        snowdepth = load_hydrologic_year(sim, year, 'Snowdepth' + var_suffix, model_outputs, glacier_ID, std_dev)
        Ptau = load_hydrologic_year(sim, year, 'Ptau' + var_suffix, model_outputs, glacier_ID, std_dev)
        SI = load_hydrologic_year(sim, year, 'Superimposedice' + var_suffix, model_outputs, glacier_ID, std_dev)
    
        for i in range(0,len(massbal)):
            massbal[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            snowmelt[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            refrozenmelt[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            netsnowmelt[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            glaciermelt[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            SImelt[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            rain[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            refrozenrain[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            rainrunoff[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            accumulation[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            snowdepth[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            Ptau[i][np.where(~np.isfinite(glacier_grid))] = np.nan
            SI[i][np.where(~np.isfinite(glacier_grid))] = np.nan

        massbal_l.append(np.nanmean(massbal,axis=(1,2)))
        snowmelt_l.append(np.nanmean(snowmelt,axis=(1,2)))
        refrozenmelt_l.append(np.nanmean(refrozenmelt,axis=(1,2)))
        netsnowmelt_l.append(np.nanmean(netsnowmelt,axis=(1,2)))
        glaciermelt_l.append(np.nanmean(glaciermelt,axis=(1,2)))
        SImelt_l.append(np.nanmean(SImelt,axis=(1,2)))
        rain_l.append(np.nanmean(rain,axis=(1,2)))
        refrozenrain_l.append(np.nanmean(refrozenrain,axis=(1,2)))
        rainrunoff_l.append(np.nanmean(rainrunoff,axis=(1,2)))
        accumulation_l.append(np.nanmean(accumulation,axis=(1,2)))
        snowdepth_l.append(np.nanmean(snowdepth,axis=(1,2)))
        Ptau_l.append(np.nanmean(Ptau,axis=(1,2)))
        SI_l.append(np.nanmean(SI,axis=(1,2)))
        gc.collect()
        
    return massbal_l, snowmelt_l, refrozenmelt_l, netsnowmelt_l, glaciermelt_l, SImelt_l, rain_l, refrozenrain_l, rainrunoff_l, accumulation_l, snowdepth_l, Ptau_l, SI_l      

def save_timeseries_as_txtfile(years,timeseries,OUTPUTPATH,modelname,variable,domain):
    ''' 
    timeseries = the list of daily timeseries for each year (e.g. massbal_kw)
    OUTPUTPATH = where the text files should be saved
    modelname (str): the name of the model (e.g. REF_MODEL, ROUNCE_DEBRIS)
    variable (str): the name of the var (e.g. massbal, totalsnowmelt)
    domain (str): what area these outputs cover (e.g. krh, allgl, kw, nonKWice)
    '''
    mb_array = np.zeros((len(years),366))
    for year in years:
        i = year-years[0]
        if len(timeseries[i]) == 366:
            mb_array[i,:] = timeseries[i]
        else:
            # if no leap year, add NaN to fill the missing day
            mb_array[i,:] = np.concatenate((timeseries[i],np.array([np.nan])))
            
    np.savetxt(os.path.join(OUTPUTPATH,modelname + '_' + variable + '_' + domain + '_' + str(years[0]) + '-' + str(years[-1]+1) + '.txt'),np.array(mb_array))

def save_timeseries(timeseries,years,output_path,variable_names,modelname,domain):
    '''
    timeseries = list of annual timeseries for each variable, produced by the fnc process_model_timeseries
    years = years to be processed (array)
    output_path = where the timeseries (.text files) will be saved
    variable_names = list of strings with variable names
    modelname = string, e.g. 'REF_MODEL'
    domain = string, describes the area the timeseries is averaged over, e.g. 'KW'
    '''
    i  = 0

    # create new directory called timeseries_1979-2022:
    new_directory = os.path.join(output_path,'timeseries_'+ str(years[0]) + '-' + str(years[-1]+1))
    os.makedirs(new_directory, exist_ok=True)
    
    for var in timeseries:
        print(variable_names[i])
        save_timeseries_as_txtfile(years,var,new_directory,modelname,variable_names[i],domain)
        i += 1
        
def process_model_distributed_fields(model_outputs, glacier_ID, sim, years, std_dev):
    '''
    model_outputs = directory with the model outputs
    glacier_ID = identifier for the glacie (string)
    sim = integer corresponding to the simulation number that is being processed
    years = years to be processed (array)
    std_dev = boolean, True if calculating standard deviation, False if just getting regular variable
    
    returns (all units = m w.e. a^-1):
    annual distributed runoff:
        net snowmelt (total snow melt minus refreezing)
        superimposed ice melt
        glacier ice melt
        rain runoff (total rain minus refrozen rain)
    
    annual distributed components:
        total snow melt
        refrozen snow melt
        total rainfall
        refrozen rainfall
        potential retention mass (Ptau)
        superimposed ice thickness
        
    annual distributed balance:
        net mass balance
        accumulation
        
    other:
        snowpack thickness
    '''
    massbal_l = []
    snowmelt_l = []
    refrozenmelt_l = []
    netsnowmelt_l = []
    glaciermelt_l = []
    SImelt_l = []
    rain_l = []
    refrozenrain_l = []
    rainrunoff_l = []
    accumulation_l = []
    snowdepth_l = []
    Ptau_l = []
    SI_l = []
    
    if std_dev == True:
        var_suffix = '_std'
    else:
        var_suffix = ''
    
    for year in years:
        print('Calculating runoff & mb variables for', year)
        # Get mass balance variable:
        massbal = load_hydrologic_year(sim, year, 'Netbalance' + var_suffix, model_outputs, glacier_ID, std_dev)
        snowmelt = load_hydrologic_year(sim, year, 'Snowmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        refrozenmelt = load_hydrologic_year(sim, year, 'Refrozenmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        netsnowmelt = load_hydrologic_year(sim, year, 'Netsnowmelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        glaciermelt = load_hydrologic_year(sim, year, 'Glaciericemelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        SImelt = load_hydrologic_year(sim, year, 'Superimposedicemelt' + var_suffix, model_outputs, glacier_ID, std_dev)
        rain = load_hydrologic_year(sim, year, 'Rain' + var_suffix, model_outputs, glacier_ID, std_dev)
        refrozenrain = load_hydrologic_year(sim, year, 'Refrozenrain' + var_suffix, model_outputs, glacier_ID, std_dev)
        rainrunoff = load_hydrologic_year(sim, year, 'Rainrunoff' + var_suffix, model_outputs, glacier_ID, std_dev)
        accumulation = load_hydrologic_year(sim, year, 'Accumulation' + var_suffix, model_outputs, glacier_ID, std_dev)
        snowdepth = load_hydrologic_year(sim, year, 'Snowdepth' + var_suffix, model_outputs, glacier_ID, std_dev)
        Ptau = load_hydrologic_year(sim, year, 'Ptau' + var_suffix, model_outputs, glacier_ID, std_dev)
        SI = load_hydrologic_year(sim, year, 'Superimposedice' + var_suffix, model_outputs, glacier_ID, std_dev)
        
        # compute sum over the hydrological year (Oct 1 - Sep 30):
        massbal_l.append(np.sum(massbal,axis=0))
        snowmelt_l.append(np.sum(snowmelt,axis=0))
        refrozenmelt_l.append(np.sum(refrozenmelt,axis=0))
        netsnowmelt_l.append(np.sum(netsnowmelt,axis=0))
        glaciermelt_l.append(np.sum(glaciermelt,axis=0))
        SImelt_l.append(np.sum(SImelt,axis=0))
        rain_l.append(np.sum(rain,axis=0))
        refrozenrain_l.append(np.sum(refrozenrain,axis=0))
        rainrunoff_l.append(np.sum(rainrunoff,axis=0))
        accumulation_l.append(np.sum(accumulation,axis=0))
        
        # for tracker variables, get value on last day of hydrological year (Sep 30):
        snowdepth_l.append(snowdepth[-1])
        Ptau_l.append(Ptau[-1])
        SI_l.append(SI[-1])
        
    return massbal_l, snowmelt_l, refrozenmelt_l, netsnowmelt_l, glaciermelt_l, SImelt_l, rain_l, refrozenrain_l, rainrunoff_l, accumulation_l, snowdepth_l, Ptau_l, SI_l

def save_distributed_fields(outputs,varnames,years,output_path,modelname,std_dev):
    i  = 0
    for output in outputs:
        print(varnames[i])
        
        # Create a directory for distributed fields:
        new_directory = os.path.join(output_path,varnames[i] + '_distributed')
        os.makedirs(new_directory, exist_ok=True)
        
        for year in years:
            print(year)
            if std_dev == False:
                np.savetxt(os.path.join(new_directory,modelname + '_' + varnames[i] + '_distributed_' + str(year) + '.txt'),output[year-years[0]])
            elif std_dev == True:
                np.savetxt(os.path.join(new_directory,modelname + '_' + varnames[i] + '_distributed_std_' + str(year) + '.txt'),output[year-years[0]])
        
        i += 1
        
# =============================================================================
# Functions for loading the formatted mbm outputs 
# =============================================================================

def load_daily_timeseries(outputs_path,modelname,varnames,domain,years):
    '''
    Load the daily timeseries for each model output saved as .text files
    (generated from the format_mbm_outputs.py script)
    '''
    
    outputs = []
    
    for var in varnames:
        print(var)
        directory = os.path.join(outputs_path,'timeseries_' + str(years[0]) + '-' + str(years[-1]+1))
        outputs.append(np.loadtxt(os.path.join(directory, modelname + '_' + var + '_' + domain + '_' + str(years[0]) + '-' + str(years[-1]+1) + '.txt')))
    
    return outputs

def load_distributed_fields(outputs_path,modelname,varnames,years):
    '''
    Load the annual distributed fields for each model output saved as .text files
    (generated from the format_mbm_outputs.py script)
    '''
    
    outputs = []
    
    for var in varnames:
        print(var)
        varlist = []
        for year in years:
            print(year)
            directory = os.path.join(outputs_path,var + '_distributed')
            varlist.append(np.loadtxt(os.path.join(directory,modelname + '_' + var + '_distributed_' + str(year) + '.txt')))
        
        outputs.append(varlist)
            
    return outputs

# =============================================================================
# Plotting functions: 
# =============================================================================

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    From https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    matplotlib.colormaps.register(cmap=newcmap)

    return newcmap

def distributed_runoff_components(years_to_plot,all_years,Sfc,runoff_component):
    '''
    Computes the runoff over a set of years and returns the distributed field
    '''
    MB_ARRAY = np.empty((len(years_to_plot),Sfc.shape[0],Sfc.shape[1]))
    for year in years_to_plot:
        #print(year)
        MB_ARRAY[year-years_to_plot[0],:,:] =  runoff_component[year-all_years[0]]

    Runoff = np.nanmean(MB_ARRAY,axis=0)
    Runoff[np.where(Runoff == 0)] = np.nan
    #print('min: ',np.nanmin(Runoff),'max: ',np.nanmax(Runoff))

    return Runoff  

def distributed_massbal_components(avg_years,all_years,Sfc,mb_component):
    MB_ARRAY = np.empty((len(avg_years),Sfc.shape[0],Sfc.shape[1]))
    for year in avg_years:
        #print(year)
        MB_ARRAY[year-avg_years[0],:,:] =  mb_component[year-all_years[0]]

    MB = np.nanmean(MB_ARRAY,axis=0)
    print('min: ',np.nanmin(MB),'max: ',np.nanmax(MB))

    return MB  


# =============================================================================
# Plotting the mass balance and runoff from the reference model: 
# =============================================================================

def distributed_mass_balance(years_to_plot,all_years,dist_massbal,contour_levels,Xgrid,Ygrid,Sfc,Catchmentoutline,KW_grid,figures_directory):
    '''
    Plots the distributed average mass balance over a set of years 
    Inputs:
        years_to_plot: array of years to be plotted (e.g. [1980, 1981, 1982] would plot
        the mass balance from Oct 1 1980-Sep 30 1982 (2 hydrological years))
        all_years: array of years that the distributed_massbal array contains
        dist_massbal = annual distributed mass balance fields
        contourf_levels = determines extent of the colorbar
        
        Xgrid, Ygrid, Sfc, Catchmentoutline, KW_grid = shapefiles defining the domain
    '''
    
    # calculate the average mass balance over the years in years_to_plot
    MB_ARRAY = np.empty((len(years_to_plot),Sfc.shape[0],Sfc.shape[1]))
    for year in years_to_plot:
        #print(year)
        MB_ARRAY[year-years_to_plot[0],:,:] = dist_massbal[year-all_years[0]]
    
    # NaN off-glacier gridcells
    mean_mb = np.nanmean(MB_ARRAY,axis=0)    
    mean_mb[Sfc==1] = np.nan
    
    kaskawulsh_mb = np.round(np.mean((np.nanmean(MB_ARRAY,axis=0))[np.isfinite(KW_grid)]),2)
    print('Kaskawulsh Glacier net mass balance =' + str(kaskawulsh_mb) + ' m w.e. a^-1 (' + str(years_to_plot[0]) + '-' + str(years_to_plot[-1]+1) + ')')  

    plt.figure(figsize=(9,5))
    plt.contourf(Xgrid,Ygrid,mean_mb,cmap='massbal',levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=np.arange(-20,5,2))
    legend.ax.set_ylabel('Mass Balance (m w.e. a$^{-1}$)', rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=16)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    plt.text(Xgrid[188,55],Ygrid[188,55],'10 km',fontsize=14,weight='bold',color='k')
    plt.hlines(Ygrid[193,49],Xgrid[193,50],Xgrid[193,100],linewidth=6,color='k')
    plt.arrow(Xgrid[162,70],Ygrid[162,70],0,3500,width=500,head_width=1800,color='k')
    plt.text(Xgrid[162,64],Ygrid[128,64],'N',fontsize=18,weight='bold',color='k')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_directory,'distributed_massbalance.pdf'),bbox_inches='tight')


def hydrograph(units,title,years_to_plot,all_years,daily_runoff_upperlim,cumu_runoff_upperlim,gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std,area_map,figures_directory):
    '''
    Plots the runoff hydrograph averaged over the years in the array years_to_plot
    
    Inputs:
        units = if units = 'mwe', y-axis units are m w.e. day-1 and m w.e. a-1, if units = 'm3', y-axis units are m3/s-1 and Gt a-1
        title = title on figure (string)
        years_to_plot: array of years to be plotted (e.g. [1980, 1981, 1982] would plot the runoff averaged from Oct 1 1980-Sep 30 1982 (2 hydrological years))
        all_years: array of years that the runoff timeseries contains
        
        daily_runoff_upperlim,cumu_runoff_upperlim = max value on left and right y-axes, respectively
        
        gl_icemelt, snowmelt_runoff, rain_runoff, superimp_icemelt,gl_icemelt_std, 
        snowmelt_runoff_std, rain_runoff_std, superimp_icemelt_std = runoff component timeseries and standard deviations
        
        area_map = grid representing the area that the timeseries was averaged over (e.g. full catchment or just the Kaskawulsh)
        
        figures_directory = path where figure will be saved
        
    Returns:
        print mean annual runoff from each component (Gt a-1) and standard deviation
        timeseries for each runoff component in m w.e. day-1
    
    '''
    
    ice_sum, snow_sum, rain_sum, SI_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    icestd_sum, snowstd_sum, rainstd_sum, SIstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in years_to_plot:
        i = year - all_years[0]
        ice_sum += np.array(gl_icemelt[i][:365])
        snow_sum += np.array(snowmelt_runoff[i][:365])
        rain_sum += np.array(rain_runoff[i][:365])
        SI_sum += np.array(superimp_icemelt[i][:365])
        
        icestd_sum += np.array(gl_icemelt_std[i][:365])
        snowstd_sum += np.array(snowmelt_runoff_std[i][:365])
        rainstd_sum += np.array(rain_runoff_std[i][:365])
        SIstd_sum += np.array(superimp_icemelt_std[i][:365])
        
    # Get mean daily values and mean daily standard deviations across all years
    ice_mean = np.array(ice_sum/len(years_to_plot))
    snow_mean = np.array(snow_sum/len(years_to_plot))
    rain_mean = np.array(rain_sum/len(years_to_plot))
    SI_mean = np.array(SI_sum/len(years_to_plot))
    total_runoff = ice_mean + snow_mean + rain_mean + SI_mean
    
    icestd_mean = np.array(icestd_sum/len(years_to_plot))
    snowstd_mean = np.array(snowstd_sum/len(years_to_plot))
    rainstd_mean = np.array(rainstd_sum/len(years_to_plot))
    SIstd_mean = np.array(SIstd_sum/len(years_to_plot))
    total_std = icestd_mean + snowstd_mean + rainstd_mean + SIstd_mean
    
    time = np.arange(0,len(total_runoff))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)
    # dc = Conversion from m w.e. / day to m3/s
    # yc = Conversion from m w.e./year to km3/year
    if units=='mwe':
        dc = 1
        yc = 1
    elif units=='m3':
        dc = area/(60*60*24) #m2/s
        yc = area/1e9 #km2

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8.5,4.73))
    ax.set_title(title,fontsize=14,loc='left',weight='bold')
    ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice runoff')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow runoff')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Refrozen ice runoff')
    ax.set_yticks(np.arange(0,451,50))
    
    plt.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='grey',alpha=0.35)
    plt.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    plt.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    plt.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    plt.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=45,fontsize=14)
    plt.xlim(0,420)

    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)

    ax0 = ax.twinx()
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    ax0.set_yticks(np.linspace(0,2.7,10))
    ax0.set_ylim(0,cumu_runoff_upperlim)
    
    if units=='mwe':
        ax.set_ylabel('Discharge (m w.e. day$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (m w.e. a$^{-1}$)',fontsize=14)
    elif units=='m3':
        ax.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
        ax0.set_ylabel('Cumulative Runoff (Gt a$^{-1}$)',fontsize=14)
    
    width=10
    plt.vlines(395, ymin=(np.cumsum(total_runoff-total_std)*yc)[-1], ymax=(np.cumsum(total_runoff+total_std)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    plt.hlines(y=(np.cumsum(total_runoff)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='k',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(ice_mean-icestd_mean)*yc)[-1], ymax=(np.cumsum(ice_mean+icestd_mean)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    plt.hlines(y=(np.cumsum(ice_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='turquoise',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(snow_mean-snowstd_mean)*yc)[-1], ymax=(np.cumsum(snow_mean+snowstd_mean)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    plt.hlines(y=(np.cumsum(snow_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='royalblue',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(rain_mean-rainstd_mean)*yc)[-1], ymax=(np.cumsum(rain_mean+rainstd_mean)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    plt.hlines(y=(np.cumsum(rain_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='deeppink',linewidth=4)
    
    plt.vlines(395, ymin=(np.cumsum(SI_mean-SIstd_mean)*yc)[-1], ymax=(np.cumsum(SI_mean+SIstd_mean)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    plt.hlines(y=(np.cumsum(SI_mean)*yc)[-1], xmin=395 - width / 2, xmax=395 + width / 2,colors='darkorange',linewidth=4)
    ax.hlines(y=-10, xmin=375 - width / 2, xmax=375 + width / 2,colors='k',linewidth=4)
    ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='Standard deviation')
    
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),loc='upper left',fontsize=14, ncol=1, borderaxespad=0.19)
    
    # Create axes for pie chart
    ax1 = plt.axes([0.05, 0.12, 0.45, 0.45])
    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']

    total_runoff_sum = np.sum(ice_mean + snow_mean + SI_mean + rain_mean)
    gl_ice_percent = round((np.sum(ice_mean)/total_runoff_sum)*100,1)
    snow_percent = round((np.sum(snow_mean)/total_runoff_sum)*100,1)
    SIice_percent = round((np.sum(SI_mean)/total_runoff_sum)*100,1)
    rain_percent = round((np.sum(rain_mean)/total_runoff_sum)*100,1)
    
    non_zero_percents = [p for p in [gl_ice_percent, snow_percent, rain_percent, SIice_percent] if p != 0]

    if all(p == 0 for p in [gl_ice_percent, snow_percent, rain_percent, SIice_percent]):
        ax1.pie([1], colors=['k'], labels=['No Data'], textprops={'fontsize': 14})
    else:
        #ax1.pie(non_zero_percents, colors=piechart_colours, textprops={'fontsize': 14}, autopct=lambda pct: '' if pct == 0 else f'{pct:.1f}%')
        ax1.pie(non_zero_percents, colors=piechart_colours)

    ax.text(151,105,str(gl_ice_percent) + '%',fontsize=14,weight='bold',color='turquoise')
    ax.text(151,80,str(snow_percent) + '%',fontsize=14,weight='bold',color='royalblue')
    ax.text(151,55,str(rain_percent) + '%',fontsize=14,weight='bold',color='deeppink')
    ax.text(151,30,str(SIice_percent) + '%',fontsize=14,weight='bold',color='darkorange')

    dates = pd.date_range(start= '2000-10-01 00:00:00',end= '2001-09-30 21:00:00',freq='1D')
    Mice_dominates = time[212:][np.where(ice_mean[212:]>snow_mean[212:])][0]
    print(Mice_dominates,str(dates[Mice_dominates])[5:10])
    #ax.text(183,400,'M$_{gl. ice}$ > M$_{snow}$\non ' + str(dates[Mice_dominates])[5:10],fontsize=14,color='darkslategrey',horizontalalignment='left')
    
    
    fig.tight_layout()
    #plt.savefig(os.path.join(figures_directory,'hydrograph.pdf'),bbox_inches='tight')

    print('glacier ice:',str(np.round((np.cumsum(ice_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('snow:',str(np.round((np.cumsum(snow_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('rain:',str(np.round((np.cumsum(rain_mean)*yc)[-1],2)) + ' Gt a$^{-1}$')
    print('refrozen ice:',str(np.round((np.cumsum(SI_mean)*yc)[-1],4)) + ' Gt a$^{-1}$')
    print('total runoff:',str(np.round(total_runoff_sum*yc,2)) + ' Gt a$^{-1}$')
    
    print('standard deviations')
    print('total runoff: ',((np.cumsum(total_runoff+total_std)*yc)[-1] - (np.cumsum(total_runoff-total_std)*yc)[-1])/2)
    print('glacier ice runoff: ',((np.cumsum(ice_mean+icestd_mean)*yc)[-1] - (np.cumsum(ice_mean-icestd_mean)*yc)[-1])/2)
    print('snow runoff: ',((np.cumsum(snow_mean+snowstd_mean)*yc)[-1] - (np.cumsum(snow_mean-snowstd_mean)*yc)[-1])/2)
    print('rain runoff: ',((np.cumsum(rain_mean+rainstd_mean)*yc)[-1] -(np.cumsum(rain_mean-rainstd_mean)*yc)[-1])/2)
    print('refrozen ice runoff: ',((np.cumsum(SI_mean+SIstd_mean)*yc)[-1] - (np.cumsum(SI_mean-SIstd_mean)*yc)[-1])/2)

    return ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent

def distributed_runoff(years_to_plot,all_years,Sfc,Xgrid,Ygrid,Catchmentoutline,gl_icemelt,netsnowmelt,rain_runoff,superimp_icemelt,figures_directory):
    '''
    Plots the distributed average runoff over a set of years 
    Inputs:
        years_to_plot: array of years to be plotted (e.g. [1980, 1981, 1982] would plot
        the mass balance from Oct 1 1980-Sep 30 1982 (2 hydrological years))
        all_years: array of years that the distributed_massbal array contains
        
        gl_icemelt,netsnowmelt,
        rain_runoff,superimp_icemelt = annual distributed runoff fields
        
        Xgrid, Ygrid, Sfc, Catchmentoutline = shapefiles defining the domain
    '''
    
    fig, axs = plt.subplots(2, 2, figsize=(9,8))
    
    # Plot 1 with inset
    icemelt = distributed_runoff_components(years_to_plot,all_years,Sfc,gl_icemelt)
    contour1 = axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18))
    axs[0, 0].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[0, 0].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[0, 0].axis('off')
    axs[0, 0].axis('equal')
    axs[0, 0].set_title('a)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider1 = make_axes_locatable(axs[0, 0])
    cax1 = divider1.append_axes("bottom", size="5%", pad=0.1)
    cb1 = plt.colorbar(contour1, cax=cax1, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb1.set_label('Runoff from glacier ice (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb1.set_ticks(np.linspace(0,20,11))
    cb1.ax.tick_params(labelsize=14)
    axs[0, 0].text(Xgrid[193,44],Ygrid[193,44],'10 km',fontsize=14,color='k')
    axs[0, 0].hlines(Ygrid[200,49],Xgrid[200,50],Xgrid[193,100],linewidth=6,color='k')
    axs[0, 0].arrow(Xgrid[170,70],Ygrid[170,70],0,3500,width=500,head_width=1800,color='k')
    axs[0, 0].text(Xgrid[135,61],Ygrid[135,61],'N',fontsize=15,weight='bold',color='k')
    
    snowmelt = distributed_runoff_components(years_to_plot,all_years,Sfc,netsnowmelt)
    contour2 = axs[0, 1].contourf(Xgrid,Ygrid,snowmelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,1,18))
    axs[0, 1].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[0, 1].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[0, 1].axis('off')
    axs[0, 1].axis('equal')
    axs[0, 1].set_title('b)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider2 = make_axes_locatable(axs[0, 1])
    cax2 = divider2.append_axes("bottom", size="5%", pad=0.1)
    cb2 = plt.colorbar(contour2, cax=cax2, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb2.set_label('Runoff from snow (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb2.set_ticks(np.linspace(0,1,6))
    cb2.ax.tick_params(labelsize=14)
    axs[0, 1].text(Xgrid[193,44],Ygrid[193,44],'10 km',fontsize=14,color='k')
    axs[0, 1].hlines(Ygrid[200,49],Xgrid[200,50],Xgrid[193,100],linewidth=6,color='k')
    axs[0, 1].arrow(Xgrid[170,70],Ygrid[170,70],0,3500,width=500,head_width=1800,color='k')
    axs[0, 1].text(Xgrid[135,61],Ygrid[135,61],'N',fontsize=15,weight='bold',color='k')
    
    rain = distributed_runoff_components(years_to_plot,all_years,Sfc,rain_runoff)
    contour3 = axs[1, 0].contourf(Xgrid,Ygrid,rain,cmap=cmocean.cm.ice_r,levels=np.linspace(0,1,28))
    axs[1, 0].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[1, 0].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[1, 0].axis('off')
    axs[1, 0].axis('equal')
    axs[1, 0].set_title('c)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider3 = make_axes_locatable(axs[1, 0])
    cax3 = divider3.append_axes("bottom", size="5%", pad=0.1)
    cb3 = plt.colorbar(contour3, cax=cax3, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb3.set_label('Runoff from rain (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb3.set_ticks(np.linspace(0,1,6))
    cb3.ax.tick_params(labelsize=14)
    axs[1, 0].text(Xgrid[193,44],Ygrid[193,44],'10 km',fontsize=14,color='k')
    axs[1, 0].hlines(Ygrid[200,49],Xgrid[200,50],Xgrid[193,100],linewidth=6,color='k')
    axs[1, 0].arrow(Xgrid[170,70],Ygrid[170,70],0,3500,width=500,head_width=1800,color='k')
    axs[1, 0].text(Xgrid[135,61],Ygrid[135,61],'N',fontsize=15,weight='bold',color='k')
    
    SI = distributed_runoff_components(years_to_plot,all_years,Sfc,superimp_icemelt)
    contour4 = axs[1, 1].contourf(Xgrid,Ygrid,SI,cmap=cmocean.cm.ice_r,levels=np.linspace(0,0.1,18))
    axs[1, 1].contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    axs[1, 1].contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    axs[1, 1].axis('off')
    axs[1, 1].axis('equal')
    axs[1, 1].set_title('d)',fontsize=14,weight='bold')
    #cb1 = plt.colorbar(axs[0, 0].contourf(Xgrid,Ygrid,icemelt,cmap=cmocean.cm.ice_r,levels=np.linspace(0,10,18)), ax=axs[0, 0])
    divider4 = make_axes_locatable(axs[1, 1])
    cax4 = divider4.append_axes("bottom", size="5%", pad=0.1)
    cb4 = plt.colorbar(contour4, cax=cax4, orientation='horizontal')
    #cb1.set_label('Value')  # Add label to colorbar
    cb4.set_label('Runoff from refrozen ice (m w.e. a$^{-1}$)', rotation=0,fontsize=14,labelpad=1)
    cb4.set_ticks(np.arange(0,0.11,0.02))
    cb4.ax.tick_params(labelsize=14)
    axs[1, 1].text(Xgrid[193,44],Ygrid[193,44],'10 km',fontsize=14,color='k')
    axs[1, 1].hlines(Ygrid[200,49],Xgrid[200,50],Xgrid[193,100],linewidth=6,color='k')
    axs[1, 1].arrow(Xgrid[170,70],Ygrid[170,70],0,3500,width=500,head_width=1800,color='k')
    axs[1, 1].text(Xgrid[135,61],Ygrid[135,61],'N',fontsize=15,weight='bold',color='k')
    plt.savefig(os.path.join(figures_directory,'distributed_runoff.pdf'),bbox_inches='tight')

# =============================================================================
# Functions to generate figures to compare the reference model with the alternative debris and accumulation models
# =============================================================================

def plot_pie_chart(ax, sizes, labels):
    piechart_colours = ['turquoise','royalblue','deeppink','darkorange']
    ax.pie(sizes, colors=piechart_colours)
    
def plot_hydrographs_template(ax,title,daily_runoff_upperlim,cumu_runoff_upperlim,percent_start,percent_gap,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, area_map):
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)

    dc = area/(60*60*24) #m2/s
    yc = area/1e9 #km2
    
    time = np.arange(0,len(total_runoff))
    ax.set_title(title,fontsize=14,loc='left',weight='bold')
    ax.plot(time,total_runoff*dc,c='k',label='Total runoff')
    ax.plot(time,ice_mean*dc,c='turquoise',label='Glacier ice runoff')    
    ax.plot(time,snow_mean*dc,c='royalblue',label='Snow runoff')
    ax.plot(time,rain_mean*dc,c='deeppink',label='Rain')
    ax.plot(time,SI_mean*dc,c='darkorange',label='Refrozen ice runoff')
    ax.set_yticks(np.arange(0,801,100))
    
    ax.fill_between(time,(total_runoff-total_std)*dc,(total_runoff+total_std)*dc,color='grey',alpha=0.35)
    ax.fill_between(time,(ice_mean-icestd_mean)*dc,(ice_mean+icestd_mean)*dc,color='turquoise',alpha=0.35)    
    ax.fill_between(time,(snow_mean-snowstd_mean)*dc,(snow_mean+snowstd_mean)*dc,color='royalblue',alpha=0.35)
    ax.fill_between(time,(rain_mean-rainstd_mean)*dc,(rain_mean+rainstd_mean)*dc,color='deeppink',alpha=0.35)
    ax.fill_between(time,(SI_mean-SIstd_mean)*dc,(SI_mean+SIstd_mean)*dc,color='darkorange',alpha=0.35)
    
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335,365])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',''],rotation=0,fontsize=14)
    ax.set_xlim(182,385)
    #ax.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
    ax.set_ylim(0,daily_runoff_upperlim)
    ax.grid()
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    ax.hlines(y=-10, xmin=375 - 10 / 2, xmax=375 + 10 / 2,colors='k',linewidth=4,label='Cumulative runoff')
    ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='Standard deviation')
    
    ax.text(185,percent_start,str(gl_ice_percent) + '%',fontsize=14,weight='bold',color='turquoise')
    ax.text(185,percent_start-percent_gap,str(snow_percent) + '%',fontsize=14,weight='bold',color='royalblue')
    ax.text(185,percent_start-(2*percent_gap),str(rain_percent) + '%',fontsize=14,weight='bold',color='deeppink')
    ax.text(185,percent_start-(3*percent_gap),str(SIice_percent) + '%',fontsize=14,weight='bold',color='darkorange')
    
    ax0 = ax.twinx()
    #if ax == ax3:
    #    pass
    #else:
    ax0.set_ylabel('Cumulative Runoff (Gt a$^{-1}$)',fontsize=14)
    width=10
    ax0.vlines(375, ymin=(np.cumsum(total_runoff-total_std)*yc)[-1], ymax=(np.cumsum(total_runoff+total_std)*yc)[-1], linewidth=width, color='k',alpha=0.3)
    ax0.hlines(y=(np.cumsum(total_runoff)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='k',linewidth=4)
    
    ax0.vlines(375, ymin=(np.cumsum(ice_mean-icestd_mean)*yc)[-1], ymax=(np.cumsum(ice_mean+icestd_mean)*yc)[-1], linewidth=width, color='turquoise',alpha=0.3)
    ax0.hlines(y=(np.cumsum(ice_mean)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='turquoise',linewidth=4)
    
    ax0.vlines(375, ymin=(np.cumsum(snow_mean-snowstd_mean)*yc)[-1], ymax=(np.cumsum(snow_mean+snowstd_mean)*yc)[-1], linewidth=width, color='royalblue',alpha=0.3)
    ax0.hlines(y=(np.cumsum(snow_mean)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='royalblue',linewidth=4)
    
    ax0.vlines(375, ymin=(np.cumsum(rain_mean-rainstd_mean)*yc)[-1], ymax=(np.cumsum(rain_mean+rainstd_mean)*yc)[-1], linewidth=width, color='deeppink',alpha=0.3)
    ax0.hlines(y=(np.cumsum(rain_mean)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='deeppink',linewidth=4)
    
    ax0.vlines(375, ymin=(np.cumsum(SI_mean-SIstd_mean)*yc)[-1], ymax=(np.cumsum(SI_mean+SIstd_mean)*yc)[-1], linewidth=width, color='darkorange',alpha=0.3)
    ax0.hlines(y=(np.cumsum(SI_mean)*yc)[-1], xmin=375 - width / 2, xmax=375 + width / 2,colors='darkorange',linewidth=4)
    
    ax0.tick_params(axis='y',labelsize=14)
    ax0.margins(x=0)
    #ax0.set_yticks(np.arange(0,2.8,0.6))
    #ax0.set_yticks(np.linspace(0,0.60,5))
    ax0.set_yticks(np.arange(0,4.6,1))
    ax0.set_ylim(0,cumu_runoff_upperlim)
    
    axins1 = inset_axes(ax,width="60%", height="60%", bbox_to_anchor=(-0.50, 0.038, 1, 1), bbox_transform=ax.transAxes)
    plot_pie_chart(axins1, sizes=[gl_ice_percent, snow_percent, rain_percent,\
 SIice_percent], labels=['A', 'B', 'C', 'D'])


def debris_sensitivity_hydrographs(years_to_plot,all_years,ref_model_runoff,debfree_model_runoff,rouncedeb_model_runoff,area_map,reference_model_outputs,debrisfree_outputs,rouncedebris_outputs):
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2,figsize=(10,6))
    
    # get runoff timeseries (units are m w.e.)
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Catchment-wide average runoff: 1980-2022 (reference model)',years_to_plot,all_years,450,2.7,ref_model_runoff[0],ref_model_runoff[1],ref_model_runoff[2],ref_model_runoff[3],ref_model_runoff[4],ref_model_runoff[5],ref_model_runoff[6],ref_model_runoff[7],area_map,reference_model_outputs) 
    plot_hydrographs_template(ax1,'a) Reference model',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, area_map)
    ax1.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
    ax2.axis('off')
    
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Debris-free model (1980-2022) ',years_to_plot,all_years,450,2.7,debfree_model_runoff[0],debfree_model_runoff[1],debfree_model_runoff[2],debfree_model_runoff[3],debfree_model_runoff[4],debfree_model_runoff[5],debfree_model_runoff[6],debfree_model_runoff[7],area_map,debrisfree_outputs)
    plot_hydrographs_template(ax3,'b) Debris-free model',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, area_map)
    ax3.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=14)
    
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Rounce et al. (2021) debris model (1980-2022)',years_to_plot,all_years,450,2.7,rouncedeb_model_runoff[0],rouncedeb_model_runoff[1],rouncedeb_model_runoff[2],rouncedeb_model_runoff[3],rouncedeb_model_runoff[4],rouncedeb_model_runoff[5],rouncedeb_model_runoff[6],rouncedeb_model_runoff[7],area_map,rouncedebris_outputs)
    plot_hydrographs_template(ax4,'c) Rounce et al. (2021)\ndebris model',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, area_map)
    
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(0.8,0.92),loc='upper right',fontsize=14, ncol=1, borderaxespad=0.19)
        
    fig.tight_layout()
    #fig.savefig(os.path.join(reference_model_outputs,'hydrographs_alt_debris_models.pdf'),bbox_inches='tight')

def ablation_debristests(years,Sfc,Xgrid,Ygrid,Zgrid,KW_grid,Catchmentoutline,abl_ref,abl_debfree,abl_rounce):
    fig = plt.figure(figsize=(8,5))
    
    gs = GridSpec(2,6, figure=fig,hspace=0,wspace=0)
    ax1 = fig.add_subplot(gs[0,:2])
    ax2 = fig.add_subplot(gs[0,2:4])
    ax3 = fig.add_subplot(gs[0,4:6])
    #ax4 = fig.add_subplot(gs[9:,0])
    ax5 = fig.add_subplot(gs[1,1:3])
    ax6 = fig.add_subplot(gs[1,3:5])

    row = 3,130
    col = 120, 329
    
    # Plot 1 with inset
    shiftedColorMap(matplotlib.cm.RdYlBu,start=0,midpoint=0.8298755186721992,stop=1,name='massbal_cmap')
    mb = distributed_massbal_components(np.arange(1980,2022),years,Sfc,abl_ref)
    mb[np.isnan(KW_grid)] = np.nan
    ax1.contourf(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],mb[row[0]:row[1],col[0]:col[1]],cmap='YlOrRd',levels=np.linspace(0,10,28))
    ax1.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Sfc[row[0]:row[1],col[0]:col[1]],levels=0,colors='k',linewidths=0.2,alpha=1)
    #ax1.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Catchmentoutline[row[0]:row[1],col[0]:col[1]],levels=1,colors='k',linewidths=0.7,alpha=1,linestyles = 'dashed')
    ax1.axis('off')
    ax1.axis('equal') 
    ax1.set_title('a)',loc='left',y=0.7,x=0.08,fontsize=10)
    
    acc = distributed_massbal_components(np.arange(1980,2022),years,Sfc,abl_debfree)
    acc[np.isnan(KW_grid)] = np.nan
    ax2.contourf(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],acc[row[0]:row[1],col[0]:col[1]],cmap='YlOrRd',levels=np.linspace(0,10,28))
    ax2.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Sfc[row[0]:row[1],col[0]:col[1]],levels=0,colors='k',linewidths=0.2,alpha=1)
    #ax2.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Catchmentoutline[row[0]:row[1],col[0]:col[1]],levels=1,colors='k',linewidths=0.7,alpha=1,linestyles = 'dashed')
    ax2.axis('off')
    ax2.axis('equal')
    ax2.set_title('b)',loc='left',y=0.7,x=0.08,fontsize=10)

    glacier = distributed_massbal_components(np.arange(1980,2022),years,Sfc,abl_rounce)
    glacier[np.isnan(KW_grid)] = np.nan
    contour3 = ax3.contourf(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],glacier[row[0]:row[1],col[0]:col[1]],cmap='YlOrRd',levels=np.linspace(0,10,28),extend='max')
    ax3.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Sfc[row[0]:row[1],col[0]:col[1]],levels=0,colors='k',linewidths=0.2,alpha=1)
    #ax3.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Catchmentoutline[row[0]:row[1],col[0]:col[1]],levels=1,colors='k',linewidths=0.7,alpha=1,linestyles = 'dashed')
    ax3.axis('off')
    ax3.axis('equal')
    ax3.set_title('c)',loc='left',y=0.7,x=0.08,fontsize=10)

    diff_nodeb = mb - acc
    shiftedColorMap(matplotlib.cm.RdBu,start=0,midpoint=0.7826086956521738,stop=1,name='massbal_diff')
    ax5.contourf(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],diff_nodeb[row[0]:row[1],col[0]:col[1]],cmap='massbal_diff',levels=np.linspace(-7.2,2,30))
    ax5.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Sfc[row[0]:row[1],col[0]:col[1]],levels=0,colors='k',linewidths=0.2,alpha=1)
    ax5.axis('off')
    ax5.axis('equal')
    ax5.set_title('d)',loc='left',y=0.7,x=0.08,fontsize=10)
    
    diff_rounce = mb - glacier
    contour6 = ax6.contourf(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],diff_rounce[row[0]:row[1],col[0]:col[1]],cmap='massbal_diff',levels=np.linspace(-7.2,2,20))
    ax6.contour(Xgrid[row[0]:row[1],col[0]:col[1]],Ygrid[row[0]:row[1],col[0]:col[1]],Sfc[row[0]:row[1],col[0]:col[1]],levels=0,colors='k',linewidths=0.2,alpha=1)
    ax6.axis('off')
    ax6.axis('equal')
    ax6.set_title('e)',loc='left',y=0.7,x=0.08,fontsize=10)   
    
    cax = fig.add_axes([0.35, 0.55, 0.33, 0.03])  # [x_position, y_position, width, height]
    cbar = fig.colorbar(contour3, cax=cax, orientation='horizontal')
    cbar.set_ticks(np.arange(0,11))
    cbar.set_label('Ablation (m w.e. a$^{-1}$)', rotation=0,fontsize=10,labelpad=1)
    cax.yaxis.tick_left()
    
    cax2 = fig.add_axes([0.35, 0.1, 0.33, 0.03])  # [x_position, y_position, width, height]
    cbar2 = fig.colorbar(contour6, cax=cax2, orientation='horizontal')
    cbar2.set_ticks(np.arange(-7,3))
    #cax2.yaxis.tick_left()
    cbar2.set_label('Ablation difference (m w.e. a$^{-1}$)', rotation=0,fontsize=10,labelpad=1)
    cax2.yaxis.tick_left()
    

    fig.tight_layout()
    #fig.savefig('D:/Model Runs/Uncertainty_Tests/debristests_distributed_ablation_JoG.pdf',bbox_inches='tight')    
    
def plot_MBcurve_template(ax,title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt, accumulation_std, refrozen_rain_std, netsnowmelt_std, superimp_icemelt_std, gl_icemelt_std, refrozen_melt, area_map):
    
    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum, refreezing_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    
    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation[i][:365])
        refrain_sum += np.array(refrozen_rain[i][:365])
        snowmelt_sum += np.array(netsnowmelt[i][:365])
        SImelt_sum += np.array(superimp_icemelt[i][:365])
        gl_melt_sum += np.array(gl_icemelt[i][:365])
        refreezing_sum += np.array(refrozen_melt[i][:365])
        
        snowfallstd_sum += np.array(accumulation_std[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std[i][:365]) 
    
    snowfall_mean = np.array(snowfall_sum/len(avg_years))
    refrain_mean = np.array(refrain_sum/len(avg_years))
    snowmelt_mean = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean = np.array(gl_melt_sum/len(avg_years))
    refreezing_mean = np.array(refreezing_sum/len(avg_years))
    
    snowmeltstd_mean = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean = np.array(gl_meltstd_sum/len(avg_years))
    
    total_accumulation = snowfall_mean + refrain_mean

    total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean
    total_ablation_std = snowmeltstd_mean + SImeltstd_mean + gl_meltstd_mean
    
    print('accumulation (m w.e. a^-1) = ',np.sum(snowfall_mean))
    print('refrozen rain (m w.e. a^-1) = ', np.sum(refrain_mean))
    print('ablation (m w.e. a^-1) = ',np.sum(total_ablation))
    print('snow melt (m w.e. a^-1) = ',np.sum(snowmelt_mean))
    print('refrozen snow melt (m w.e. a^-1) = ',np.sum(refreezing_mean))
    
    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9 # Gt
    area = 1 # mw.e.
    
    ax.set_title(title,fontsize=11,loc='left',y=0.87,x=0.02)
    ax.plot(time,total_accumulation,c='mediumblue',label='Accumulation')    
    ax.plot(time,-total_ablation,c='red',label='Ablation')    
    #ax.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    ax.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','','Dec','','Feb','','Apr','','Jun','','Aug',''],rotation=0,fontsize=10)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    #ax.set_ylim(-0.03,0.01)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=10)
    ax.margins(x=0)
    #ax.set_xlim(182,365)
    
    massbal = total_accumulation - total_ablation

    transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    
    if title[0:2] == 'c)':
        ax0.set_ylabel('Cumulative Balance (m w.e. a$^{-1}$)',fontsize=10)
        ax0.set_yticks(np.arange(-0.9,0.91,0.3))
        ax.set_yticklabels([])
        ax0.plot(time,np.cumsum(massbal)*area,c='k',label='$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    elif title[0:2] == 'b)':
        ax0.set_yticks(np.arange(-0.9,0.91,0.3))
        ax0.set_yticklabels([])
        ax.set_yticklabels([])
        #ax.legend(fontsize=10,loc='upper right')
        ax.plot(time,np.cumsum(massbal)*area-100,c='k',label='Cumulative balance',linewidth=3)
        ax0.plot(time,np.cumsum(massbal)*area,c='k',label='$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    elif title[0:2] == 'a)':
        ax0.set_yticks(np.arange(-0.9,0.91,0.3))
        ax0.set_yticklabels([])
        ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=10)
        ax0.plot(time,np.cumsum(massbal)*area,c='k',label='$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    
    ax0.legend(fontsize=11,loc='lower left')
    
    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    #ax0.set_ylim(-2.7,0.9)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=10)
    
    
def accumulation_sensitivity_tests(years_to_plot,all_years,ref_model_runoff,uncorracc_model_runoff,awsbias_model_runoff,ref_model_mb,uncorracc_mb,awsbias_mb,runoff_area_map,mb_area_map,reference_model_outputs,uncorracc_model_outputs,awsbias_model_outputs):
    fig = plt.figure(figsize=(11, 6))
    #fig = plt.figure(figsize=(11, 5.625))
    
    # Create a GridSpec object with 2 rows and 3 columns
    gs = GridSpec(2, 3, figure=fig,height_ratios=[0.7,1]) # top row looks good at 3.125 inches tall
    
    # get runoff timeseries (units are m w.e.)
    ax1 = fig.add_subplot(gs[1, 0])
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Catchment-wide average runoff: 1980-2022 (reference model)',years_to_plot,all_years,450,2.7,ref_model_runoff[0],ref_model_runoff[1],ref_model_runoff[2],ref_model_runoff[3],ref_model_runoff[4],ref_model_runoff[5],ref_model_runoff[6],ref_model_runoff[7],runoff_area_map,reference_model_outputs) 
    plot_hydrographs_template(ax1,'d)',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, runoff_area_map)
    ax1.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=10)
    
    ax2 = fig.add_subplot(gs[1, 1])
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Catchment-wide average runoff: 1980-2022 (reference model)',years_to_plot,all_years,450,2.7,uncorracc_model_runoff[0],uncorracc_model_runoff[1],uncorracc_model_runoff[2],uncorracc_model_runoff[3],uncorracc_model_runoff[4],uncorracc_model_runoff[5],uncorracc_model_runoff[6],uncorracc_model_runoff[7],runoff_area_map,uncorracc_model_outputs) 
    plot_hydrographs_template(ax2,'e)',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, runoff_area_map)
    #ax2.set_ylabel('Discharge (m$^3$ s$^{-1}$)',fontsize=10)
    #ax2.set_yticks(np.arange(0,451,100))
    
    ax3 = fig.add_subplot(gs[1, 2])
    ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, \
    rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent,\
     SIice_percent = hydrograph('m3','Catchment-wide average runoff: 1980-2022 (reference model)',years_to_plot,all_years,450,2.7,awsbias_model_runoff[0],awsbias_model_runoff[1],awsbias_model_runoff[2],awsbias_model_runoff[3],awsbias_model_runoff[4],awsbias_model_runoff[5],awsbias_model_runoff[6],awsbias_model_runoff[7],runoff_area_map,awsbias_model_outputs) 
    plot_hydrographs_template(ax3,'f)',450,2.7,160,40,ice_mean,snow_mean, rain_mean, SI_mean, total_runoff, icestd_mean, snowstd_mean, rainstd_mean, SIstd_mean, total_std, gl_ice_percent, snow_percent, rain_percent, SIice_percent, runoff_area_map)
    #ax3.set_yticks(np.arange(0,451,100))
    
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label.values(), by_label.keys(),loc='upper right',fontsize=10, ncol=3, borderaxespad=0.19,bbox_to_anchor=(0.92,1.01))
    
    ax4 = fig.add_subplot(gs[0, 0])
    plot_MBcurve_template(ax4,'a)',np.arange(1980,2022),all_years,0.03,0.9, ref_model_mb[0], ref_model_mb[1], ref_model_mb[2], ref_model_mb[3], ref_model_mb[4], ref_model_mb[5], ref_model_mb[6], ref_model_mb[7], ref_model_mb[8], ref_model_mb[9],ref_model_mb[10],mb_area_map)
    
    ax5 = fig.add_subplot(gs[0, 1])
    plot_MBcurve_template(ax5,'b)',np.arange(1980,2022),all_years,0.03,0.9, uncorracc_mb[0], uncorracc_mb[1], uncorracc_mb[2], uncorracc_mb[3], uncorracc_mb[4], uncorracc_mb[5], uncorracc_mb[6], uncorracc_mb[7], uncorracc_mb[8], uncorracc_mb[9],uncorracc_mb[10],mb_area_map)
    
    ax6 = fig.add_subplot(gs[0, 2])
    plot_MBcurve_template(ax6,'c)',np.arange(1980,2022),all_years,0.03,0.9, awsbias_mb[0], awsbias_mb[1], awsbias_mb[2], awsbias_mb[3], awsbias_mb[4], awsbias_mb[5], awsbias_mb[6], awsbias_mb[7], awsbias_mb[8], awsbias_mb[9],awsbias_mb[10],mb_area_map)
    
    handles, labels = ax5.get_legend_handles_labels()
    by_label2 = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    fig.legend(by_label2.values(), by_label2.keys(),loc='upper right',fontsize=10, ncol=2, borderaxespad=0.19,bbox_to_anchor=(0.4,1.01))
    
    
    fig.tight_layout()
    #fig.savefig('D:/Model Runs/Uncertainty_Tests/AccumulationTests_Hydrographs_Catchmentwide_1980-2022_JoG.pdf',bbox_inches='tight')

def massbalance_timeseries_average(title,avg_years,all_years,daily_mb_abs_lim,cumu_mb_abs_lim, accumulation, refrozen_rain, netsnowmelt, superimp_icemelt, gl_icemelt, accumulation_std, refrozen_rain_std, netsnowmelt_std, superimp_icemelt_std, gl_icemelt_std, refrozen_melt, area_map):
    
    snowfall_sum, refrain_sum, snowmelt_sum, SImelt_sum, gl_melt_sum, refreezing_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))
    snowfallstd_sum, refrainstd_sum, snowmeltstd_sum, SImeltstd_sum, gl_meltstd_sum = np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365)), np.zeros((365))

    for year in avg_years:
        i = year - all_years[0]
        snowfall_sum += np.array(accumulation[i][:365])
        refrain_sum += np.array(refrozen_rain[i][:365])
        snowmelt_sum += np.array(netsnowmelt[i][:365])
        SImelt_sum += np.array(superimp_icemelt[i][:365])
        gl_melt_sum += np.array(gl_icemelt[i][:365])
        refreezing_sum += np.array(refrozen_melt[i][:365])
        
        snowfallstd_sum += np.array(accumulation_std[i][:365])
        refrainstd_sum += np.array(refrozen_rain_std[i][:365])
        snowmeltstd_sum += np.array(netsnowmelt_std[i][:365])
        SImeltstd_sum += np.array(superimp_icemelt_std[i][:365])
        gl_meltstd_sum += np.array(gl_icemelt_std[i][:365]) 
        
    snowfall_mean = np.array(snowfall_sum/len(avg_years))
    refrain_mean = np.array(refrain_sum/len(avg_years))
    snowmelt_mean = np.array(snowmelt_sum/len(avg_years))
    SImelt_mean = np.array(SImelt_sum/len(avg_years))
    gl_melt_mean = np.array(gl_melt_sum/len(avg_years))
    refreezing_mean = np.array(refreezing_sum/len(avg_years))
      
    snowfallstd_mean = np.array(snowfallstd_sum/len(avg_years))
    refrainstd_mean = np.array(refrainstd_sum/len(avg_years))
    snowmeltstd_mean = np.array(snowmeltstd_sum/len(avg_years))
    SImeltstd_mean = np.array(SImeltstd_sum/len(avg_years))
    gl_meltstd_mean = np.array(gl_meltstd_sum/len(avg_years))

    total_accumulation = snowfall_mean + refrain_mean
    total_accumulation_std = snowfallstd_mean + refrainstd_mean
    
    total_ablation = snowmelt_mean + SImelt_mean + gl_melt_mean
    total_ablation_std = snowmeltstd_mean + SImeltstd_mean + gl_meltstd_mean
    
    print('accumulation (m w.e. a^-1) = ',np.sum(snowfall_mean))
    print('refrozen rain (m w.e. a^-1) = ', np.sum(refrain_mean))
    print('ablation (m w.e. a^-1) = ',np.sum(total_ablation))
    print('snow melt (m w.e. a^-1) = ',np.sum(snowmelt_mean))
    print('refrozen snow melt (m w.e. a^-1) = ',np.sum(refreezing_mean))

    dates = pd.date_range(start= str(2008) + '-10-01 00:00:00',end= str(2008+1) + '-09-30 21:00:00',freq='1D')   
    time = np.arange(0,len(dates))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9

    # Plot the average
    #fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,3.5))
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(7,4))
    #ax.set_title(title + str(avg_years[0])+'-'+str(avg_years[-1]+1),fontsize=14)
    ax.set_ylabel('Mass Change (m w.e. d$^{-1}$)',fontsize=15)
    ax.plot(time,total_accumulation,c='mediumblue',label='Accumulation')    
    ax.plot(time,-total_ablation,c='red',label='Ablation')    
    plt.fill_between(time,(total_accumulation-total_accumulation_std),(total_accumulation+total_accumulation_std),color='mediumblue',alpha=0.35)    
    plt.fill_between(time,(-total_ablation-total_ablation_std),(-total_ablation+total_ablation_std),color='red',alpha=0.35)
    ax.set_xticks(ticks=[0,31,61,92,123,151,182,212,243,273,304,335])
    ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'],rotation=45,fontsize=14)
    ax.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax.grid()
    ax.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax.tick_params(axis='y',labelsize=14)
    ax.margins(x=0)
    
    massbal = total_accumulation - total_ablation
    max_mb = (total_accumulation + total_accumulation_std) - (total_ablation - total_ablation_std)
    min_mb = (total_accumulation - total_accumulation_std) - (total_ablation + total_ablation_std)

    transition_date = dates[50:][np.where(np.cumsum(massbal)[50:] <=0)[0][0]]
    ax0 = ax.twinx()
    ax0.plot(time,np.cumsum(massbal)*area,c='k',label='Cumulative balance\n$\dot{B}$ = 0 on ' + str(transition_date)[5:10],linewidth=3)
    print('mass bal (Gt yr-1) =',(np.cumsum(massbal)*area)[-1])
    print('std dev (Gt yr-1) =',((np.cumsum(max_mb)*area)[-1]) - ((np.cumsum(min_mb)*area)[-1]))
    #plt.fill_between(time,np.cumsum(min_mb),np.cumsum(max_mb),color='k',alpha=0.35)    
    #ax.fill_between(time,-12,-10,color='grey',alpha=0.35,label='std. dev.')

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass Change (Gt yr$^{-1}$)',fontsize=15)
    ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=15,loc='lower left')
        
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax.legend(by_label.values(), by_label.keys(),fontsize=15,loc='upper left', ncol=1, borderaxespad=0.19)
    #fig.patch.set_facecolor('#f8f5f0')
    fig.tight_layout()
    print('cumulative massbal (m w.e. a) = ',(np.cumsum(massbal))[-1])
    print('max winter balance (m w.e. a) = ',np.max(np.cumsum(massbal)))
    print('date of max winter balance = ',dates[np.where(np.cumsum(massbal) == np.max(np.cumsum(massbal)))])
    #fig.savefig('D:\Model Runs\REF_MODEL\Plots\AnnualMB_1998-1999.pdf',bbox_inches='tight')

def cumulative_and_annual_massbal(title,years,plotting_years,daily_mb_abs_lim,cumu_mb_abs_lim,mb_abs_lim,accumulation,refrozen_rain,gl_icemelt,netsnowmelt,SImelt,area_map):
    year_i_index = plotting_years[0] - years[0]
    year_f_index = len(years) - (years[-1] - plotting_years[-1])
    
    all_accumulation = np.concatenate(accumulation[year_i_index:year_f_index], axis=0) + np.concatenate(refrozen_rain[year_i_index:year_f_index], axis=0)
    all_ablation = np.concatenate(gl_icemelt[year_i_index:year_f_index], axis=0) + np.concatenate(netsnowmelt[year_i_index:year_f_index], axis=0) + np.concatenate(SImelt[year_i_index:year_f_index], axis=0)
    
    acc_final = all_accumulation[~np.isnan(all_accumulation)]
    abl_final = all_ablation[~np.isnan(all_ablation)]
    
    massbal = acc_final - abl_final    
    
    dates = pd.date_range(start= str(plotting_years[0]) + '-10-01 00:00:00',end= str(plotting_years[-1]+1) + '-09-30 21:00:00',freq='1D')   
    filler_dates = len(dates) - (((plotting_years.repeat(365).reshape(len(plotting_years),365)) + np.linspace(0,0.995,365)).flatten()).shape[0]
    datess = np.concatenate((((plotting_years.repeat(365).reshape(len(plotting_years),365)) + np.linspace(0,0.995,365)).flatten(),np.arange(plotting_years[-1]+1,plotting_years[-1]+1+(filler_dates*0.002733),0.002733)))
    
    area = np.where(np.isfinite(area_map))[0].shape[0]*(200*200)/1e9
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1,figsize=(10,7))
    #ax1.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax1.set_ylabel('Mass Change\n(m w.e. d$^{-1}$)',fontsize=14)
    ax1.plot(datess,acc_final,c='dodgerblue',label='Accumulation')    
    ax1.plot(datess,-abl_final,c='crimson',label='Ablation')  
    ax1.set_ylim(-daily_mb_abs_lim,daily_mb_abs_lim)
    ax1.grid()
    ax1.axhline(y=0,xmin=0,xmax=len(dates),linestyle='--',c='k')
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)
    #ax1.margins(x=0)
    ax1.text(1980,0.06,'a)',fontsize=14,weight='bold')
    ax1.set_xlim(1979,2022)
    
    ax0 = ax1.twinx()
    ax0.plot(datess,np.cumsum(massbal)*area,c='k',linewidth=3,label='Cumulative Mass Balance')
    print((np.cumsum(massbal)*area)[-1])

    ax0.set_ylim(-cumu_mb_abs_lim,cumu_mb_abs_lim)
    ax0.set_ylabel('Cumulative Mass\nBalance (Gt)',fontsize=14)
    #ax0.margins(x=0)
    ax0.tick_params(axis='y',labelsize=14)
    ax0.legend(fontsize=14,loc='lower left')
        
    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    
    ax1.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper right', ncol=2, borderaxespad=0.19)
    
    all_accumulation = np.nansum(accumulation[year_i_index:year_f_index], axis=1) + np.nansum(refrozen_rain[year_i_index:year_f_index], axis=1)
    all_ablation = np.nansum(gl_icemelt[year_i_index:year_f_index], axis=1) + np.nansum(netsnowmelt[year_i_index:year_f_index], axis=1) + np.nansum(SImelt[year_i_index:year_f_index], axis=1)
    massbal = all_accumulation - all_ablation 
    
    #ax2.set_title(title + str(plotting_years[0])+'-'+str(plotting_years[-1]+1),fontsize=14)
    ax2.set_ylabel('Mass Balance\n(m w.e. a$^{-1}$)',fontsize=14)
    ax2.bar(plotting_years,all_accumulation,color='dodgerblue',label='Accumulation',zorder=10)    
    ax2.bar(plotting_years,-all_ablation,color='crimson',label='Ablation',zorder=10)  
    ax2.set_ylim(-mb_abs_lim,mb_abs_lim)
    ax2.grid(zorder=-10)
    ax2.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')
    ax2.tick_params(axis='y',labelsize=14)
    ax2.tick_params(axis='x',labelsize=14)
    #ax2.margins(x=0.01)
    ax2.text(1980,2,'b)',fontsize=14,weight='bold')
    ax2.set_xlim(1979,2022)
    
    ax3 = ax2.twinx()
    #ax3.bar(plotting_years,massbal*area*0,color='white',linewidth=3,label='Mass Balance',zorder=-20)
    ax3.set_ylim(-mb_abs_lim*area,mb_abs_lim*area)
    ax3.set_yticks(np.round([-1*area,-2*area,0,1*area,2*area],1))
    ax3.tick_params(axis='y',labelsize=14)
    ax3.set_ylabel('Mass Balance\n(Gt a$^{-1}$)',fontsize=14)
    
    ax2.bar(plotting_years,massbal,color='k',linewidth=3,label='Mass Balance',zorder=30)
    #ax3.axhline(y=0,xmin=0,xmax=len(plotting_years),linestyle='--',c='k')    
    
    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #plt.legend(by_label.values(), by_label.keys())
    ax2.legend(by_label.values(), by_label.keys(),fontsize=14,loc='upper right', ncol=3, borderaxespad=0.19)
    fig.tight_layout()
    
    return all_accumulation, -all_ablation, massbal

def distributed_var(plotting_years,years,var,contour_levels,label,cmap,Xgrid,Ygrid,Sfc,Catchmentoutline):
    i_start = plotting_years[0] - years[0]
    i_end = plotting_years[-1] - years[0] + 1
    
    plt.figure(figsize=(9,5))
    plt.title(str(plotting_years[0]) + '-' + str(plotting_years[-1]+1),fontsize=12)
    #plt.title(model_name,fontsize=14)
    plt.contourf(Xgrid,Ygrid,np.nanmean(np.array(var[i_start:i_end]),axis=0),cmap=cmap,levels=contour_levels)
    plt.axis('equal')
    legend = plt.colorbar(ticks=contour_levels)
    legend.ax.set_ylabel(label, rotation=270,fontsize=14,labelpad=25)
    plt.xlabel('Easting (m)',fontsize=14)
    plt.ylabel('Northing (m)',fontsize=14)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    legend.ax.tick_params(labelsize=16)
    plt.contour(Xgrid,Ygrid,Sfc,levels=0,colors='k',linewidths=0.5,alpha=0.8)
    #plt.contour(Xgrid,Ygrid,np.nanmean(MB_ARRAY,axis=0),levels=0,colors='k',linewidths=1,alpha=1)
    plt.contour(Xgrid,Ygrid,Catchmentoutline,levels=1,colors='k',linewidths=0.9,alpha=1,linestyles = 'dashed')
    #plt.text(567000,6710000,str(avg_years[0]) + '--' + str(avg_years[-1]+1) + '\nmass balance\n= ' + str(kaskawulsh_mb) + ' m w.e. a$^{-1}$',fontsize=15)    
    plt.text(Xgrid[188,55],Ygrid[188,55],'10 km',fontsize=14,weight='bold',color='k')
    plt.hlines(Ygrid[193,49],Xgrid[193,50],Xgrid[193,100],linewidth=6,color='k')
    plt.arrow(Xgrid[162,70],Ygrid[162,70],0,3500,width=500,head_width=1800,color='k')
    plt.text(Xgrid[162,64],Ygrid[128,64],'N',fontsize=18,weight='bold',color='k')
    # fig.patch.set_facecolor('#f8f5f0')
    plt.axis('off')
    plt.tight_layout()