# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:40:57 2023

Note: the time-weighted averaging is specific to this area (based on the synoptic time scale in the St. Elias Mountains)
and is specific to the set of satellite images used here. This code is cannot be used to calculate snowline scores 
or perform stage II of tuning for other glaciers without modifying the time-weighted averaging. 

This script also plots figures 5.3, 5.4, 5.15 from KR's thesis

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import heapq
from scipy.stats import norm

model_name = 'REF_MODEL' # string assigned to this run of the model to help keep track of results

# parameters and modelled mass balances from sims that passed stage I:
stageI_params_file = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageI/tuning_stageI_params_REF_MODEL.csv' # Path to csv file containing aice, asnow, MF used for tuning
modelled_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageI/tuning_stageI_mb_REF_MODEL.csv' #  csv file containing the net mb for each sim (calculated using calulate_mb.py)

# observed/modelled snowline scores:
observed_snowline_rasters = '/home/krobin/projects/def-gflowers/krobin/MassBalanceModel/Inputs/KRH/Snowlines' # directory with rasterized observed snowlines (.npy files)
modelled_snowline_scores = '/home/krobin/scratch/ModelRuns/KRH/VAAtest2'

# results of tuning stage II:
output_stageII_params = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageII' 

# optional: compare reference model with tuned parameters from the alternative debris and accumulation tests (Robinson et al. 2024)
compare_alt_model_params = True
rounce_debris_finalparams = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageII/tuning_stageII_params_ROUNCE_DEBRIS.csv'
debris_free_finalparams = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageII/tuning_stageII_params_NO_DEBRIS.csv'
uncorrected_acc_finalparams = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageII/tuning_stageII_params_UNCORRECTED_ACC.csv'
aws_bc_finalparams = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/stageII/tuning_stageII_params_AWS_BIAS.csv'

# =============================================================================
# Load the tuning results from stage I:
# =============================================================================

params = np.loadtxt(stageI_params_file,skiprows=1,delimiter=',') 
aice = params[:,1]
asnow = params[:,2]
MF = params[:,3]

stageI_mb = np.loadtxt(modelled_mb,skiprows=1,delimiter=',')
mb = stageI_mb[:,2]

simID = np.arange(0,len(mb)) # ordered from 0 to len(mb)
stageIresults_simID = stageI_mb[:,1] # original sim IDs from the initial 10k parameters

# =============================================================================
# Remove sims where aice < asnow
# =============================================================================

mb_aice_geq_asnow = np.delete(mb,np.where(aice<asnow))
simID_aice_geq_asnow = np.delete(simID,np.where(aice<asnow))

target_mb = -0.46
tagret_mb_stddev = 0.17

# =============================================================================
# Get simIDs and param values for all sims (excluding ones with aice <asnow) where the 2007-18 MB is within -0.46 +- 3std.dev
# =============================================================================

simID_passing = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(3*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(3*tagret_mb_stddev))))]
mb_passing = mb_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (target_mb-(3*tagret_mb_stddev))) & (mb_aice_geq_asnow <= (target_mb+(3*tagret_mb_stddev))))]
aice_passing = aice[simID_passing]
asnow_passing = asnow[simID_passing]
MF_passing = MF[simID_passing]

stageIresults_simID_passing = stageIresults_simID[simID_passing]

# =============================================================================
# Get list of dates where there is an observed snowline raster:
# =============================================================================

all_snow_depth_dates = []
for file in os.listdir(observed_snowline_rasters):
    if file.endswith('.npy'):
        date = pd.Timestamp(year=int(file[15:19]), month=int(file[20:22]), day=int(file[23:25]))
        all_snow_depth_dates.append(date)

# =============================================================================
# Calculate the time-averaged snowline score for each simulation
# =============================================================================

all_sim_scores = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_scores[:] = np.nan

all_sim_weights = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_weights[:] = np.nan

all_sim_SA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_SA[:] = np.nan

all_sim_SW = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_SW[:] = np.nan

all_sim_CA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_CA[:] = np.nan

all_sim_NA = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_NA[:] = np.nan

all_sim_TR = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_TR[:] = np.nan

all_sim_weightedscores = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_weightedscores[:] = np.nan

all_sim_norm = np.zeros((len(mb_passing),len(all_snow_depth_dates)))
all_sim_norm[:] = np.nan

sim = 0
for i in simID_passing:
    print(i)
    scores = np.loadtxt(os.path.join(modelled_snowline_scores,'sim' + str(sim) + '_snowlinescores.txt'))
    #scores = np.loadtxt(os.path.join(modelled_snowline_scores,'StageII_scores_sim' + str(sim) + '.txt'))
    all_sim_scores[sim,:] = scores[0]
    all_sim_weights[sim,:] = scores[1]
    all_sim_SA[sim,:] = scores[2]
    all_sim_SW[sim,:] = scores[3]
    all_sim_CA[sim,:] = scores[4]
    all_sim_NA[sim,:] = scores[5]
    all_sim_TR[sim,:] = scores[6]
    all_sim_weightedscores[sim,:] = (scores[0]*scores[1])
    all_sim_norm[sim,:] = scores[7]
    sim += 1
    
# Implement the time-dependent averaging (this is specific to the Kaskawulsh and this selection of satellite images):
def get_final_timeaveraged_scores(score):
    
    time_averaged_score = np.array((np.arange(0,len(mb_passing))))
    
    for day in all_snow_depth_dates:
        #print(day)
        if day == pd.Timestamp(2016,7,4):
            pass
        elif day == pd.Timestamp(2017,7,24):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,7,24))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,7,25))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,7,25):
            pass
        
        elif day == pd.Timestamp(2017,8,6):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,8,6))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,8,8))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,8,8):
            pass
        
        elif day == pd.Timestamp(2017,9,13):
            pass
        
        elif day == pd.Timestamp(2017,9,20):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,9,20))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2017,9,23))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2017,9,23):
            pass  
        
        elif day == pd.Timestamp(2018,6,10):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,10))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,14))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,6,14):
            pass  
        
        elif day == pd.Timestamp(2018,6,25):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,25))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,6,27))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,6,27):
            pass  
        
        elif day == pd.Timestamp(2018,7,22):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,22))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,24))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,25))[0][0]])/3
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,7,24):
            pass  
        elif day == pd.Timestamp(2018,7,25):
            pass 
        
        elif day == pd.Timestamp(2018,7,30):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,7,30))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,1))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,8,1):
            pass 
        
        elif day == pd.Timestamp(2018,8,18):
            average_score = (score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,18))[0][0]] + score[:,np.where(np.array(all_snow_depth_dates)==pd.Timestamp(2018,8,23))[0][0]])/2
            time_averaged_score = np.column_stack((time_averaged_score,average_score))
        elif day == pd.Timestamp(2018,8,23):
            pass  
        
        else:
            time_averaged_score = np.column_stack((time_averaged_score,score[:,np.where(np.array(all_snow_depth_dates)==day)][:,0,0]))
        
    return time_averaged_score
        
time_avg_score = get_final_timeaveraged_scores(all_sim_weightedscores)[:,1:]
final_score = np.mean(time_avg_score,axis=1)
 
# Calculate the maximum possible score if all images were perfectly matches by the model (all image scores = 1)
max_possible_score = get_final_timeaveraged_scores(np.ones((all_sim_scores.shape))*all_sim_weights)[:,1:]
maximum_final_score = np.max(np.mean(max_possible_score,axis=1))
       
#normalize the final scores:
final_scores_norm = np.array(final_score)/maximum_final_score

# Get final scores for each tributary:
time_avg_score = get_final_timeaveraged_scores(all_sim_SA)[:,1:]
final_score_SA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_SW)[:,1:]
final_score_SW = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_CA)[:,1:]
final_score_CA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_NA)[:,1:]
final_score_NA = np.nanmean(time_avg_score,axis=1)

time_avg_score = get_final_timeaveraged_scores(all_sim_TR)[:,1:]
final_score_TR = np.nanmean(time_avg_score,axis=1)


# =============================================================================
# Select the 100 best snowline scores from a normal distribution imposed on the binned modelled mass balances
# =============================================================================

N = int(np.round(np.sqrt(len(mb_passing))))
mb_normaldist_bins =  np.linspace(-0.46-(3*0.17),-0.46+(3*0.17),N)
bin_centers = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3.1*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0]))

#simID_stageI = simID_aice_geq_asnow[np.where((mb_aice_geq_asnow >= (-0.46-(3*0.17))) & (mb_aice_geq_asnow <= (-0.46+(3*0.17))))]
simID_stageI = np.arange(0,len(mb_passing))

simID_under_normaldist = []
mb_under_normaldist = []
snowlines_under_normaldist = []
N_sims_per_bin = []
for bin_edge in range(0,len(mb_normaldist_bins)-1):
    # Define the normal distribution:
    #x = np.arange(-0.46-(3*0.17)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,-0.46+(3*0.17)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0])) # aligns with the center of each bin
    p = norm.pdf(bin_centers, -0.46,0.17)
    p_scaled =  p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][4]
    p_scaled =  p/np.max(p)*9.6
    
    # Define edges of the bin
    left_edge = np.round(mb_normaldist_bins[bin_edge],2) # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = np.round(mb_normaldist_bins[bin_edge+1],2)
    
    left_edge = mb_normaldist_bins[bin_edge] # bins are defined by [left_edge,right_edge), except the last bin where the right_edge is also included
    right_edge = mb_normaldist_bins[bin_edge+1]
    #print(bin_edge,left_edge,right_edge)
    
    bin_center = np.round((left_edge + right_edge)/2,2)
    print(bin_edge,bin_center)
    
    # Get number of sims in this bin allowed by the normal distribution
    N = int(np.round((p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]])))
    print(bin_center,p_scaled[np.where(np.round(bin_centers,2)==bin_center)[0][0]],N)
    N_sims_per_bin.append(N)
    
    if N == 0:
        pass
    else:
        
        # Get the sim IDs and snowline scores corresponding to the mass balance values within each bin 
        simID_bin = simID_stageI[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        snowline_scores_bin = final_scores_norm[np.where((mb_passing >= (left_edge)) & (mb_passing < (right_edge)))]
        
        # Sample N unique elements from the bin
        if N > len(simID_bin):
            N = len(simID_bin)
            
        # Instead of randomly sampling from inside the bin, get the N best snowline scores from inside the bin
        best_scores = (heapq.nlargest(N, snowline_scores_bin))
        # get Sim IDs corresponding to the best scores
        sims = simID_bin[np.where(snowline_scores_bin >= np.min(best_scores))]
        
        # save the sample to the normal distirbution lists
        simID_under_normaldist.extend(sims)
        mb_under_normaldist.extend(list(mb_passing[sims]))
        snowlines_under_normaldist.extend(list(final_scores_norm[sims]))
  
# Save final sims:
aice_final = aice_passing[simID_under_normaldist]
asnow_final = asnow_passing[simID_under_normaldist]
MF_final = MF_passing[simID_under_normaldist]
simID_final = simID_passing[simID_under_normaldist]
stageIresults_simID_final = stageIresults_simID_passing[simID_under_normaldist]

# =============================================================================
# Save the final params for this model
# =============================================================================
  
d = {'aice': list(aice_final), 'asnow': list(asnow_final), 'MF': list(MF_final)}
df = pd.DataFrame(data=d)
df.to_csv(os.path.join(output_stageII_params,'tuning_stageII_params_' + model_name + '.csv'))

d2 = {'simID': list(stageIresults_simID_final), 'MB': list(mb_under_normaldist), 'Snowlinescore': list(snowlines_under_normaldist)}
df2 = pd.DataFrame(data=d)
df2.to_csv(os.path.join(output_stageII_params,'tuning_stageII_mb_' + model_name + '.csv'))

# =============================================================================
# Plot the final parameter selection from the binned modelled mass balances (fig 5.3 from KR thesis)
# =============================================================================

fig = plt.figure(figsize=(7.5,3.5))
plt.grid(zorder=20)
#plt.ylim(0,170)
plt.text(-1,67,'Site-specific reference model',weight='bold',fontsize=12)
plt.hist(mb_passing,mb_normaldist_bins,rwidth=0.9,color='slategrey',zorder=10,label='Simulations within $\pm$ 3$\sigma$ of $\dot{B}_{obs}$,\n(Total = ' + str(len(mb_passing)) +')')
plt.hist(mb_under_normaldist,mb_normaldist_bins,rwidth=0.9,color='navy',zorder=20,label='Final simulations, (Total = 100)')
plt.ylabel('Frequency',fontsize=12)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(bin_centers,np.round(bin_centers,2),fontsize=12,rotation=45)
x = np.arange(-0.98,0.11,0.02) # aligns with the center of each bin
p = norm.pdf(x, -0.46,0.17)
plt.ylim(0,65)
#plt.plot(x, p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][4], 'k', linewidth=4,zorder=30)
plt.plot(x, p/np.max(p)*9.6, 'k', linewidth=4,zorder=30,label='Target distribution')
fig.legend(fontsize=12,ncol=1,bbox_to_anchor=(0.98,0.88))
plt.tight_layout()
plt.savefig(os.path.join(output_stageII_params,'Param_selection_from_MB_distribution.pdf'),bbox_inches='tight')

# =============================================================================
# Scatter plot of modelled mass bal vs snowline score (fig 5.4 from KR thesis)
# =============================================================================

fig = plt.figure(figsize=(4.8,4))
plt.text(0.80,0.13,'Site-specific reference model',weight='bold',fontsize=14)
#plt.title('Final param selection\nTotal = ' + str(len(mb_under_normaldist)),fontsize=14)
plt.scatter(final_scores_norm,mb_passing,c='slategrey',s=20,label='Simulations within\n$\pm$ 3$\sigma$ of $\dot{B}_{obs}$,\n(Total = ' + str(len(mb_passing))+ ')')
plt.scatter(snowlines_under_normaldist,mb_under_normaldist,c='navy',edgecolor='navy',s=20,label='Final simulations,\n(Total = 100)')
#plt.scatter(final_scores_norm[np.where(aice_passing>asnow_passing)],mb_passing[np.where(aice_passing>asnow_passing)],c='darkblue',label='a$_{ice}$ $\geq$ a$_{snow}$')
plt.xlabel('Snowline Score (normalized) (-)',fontsize=14)
plt.ylabel('2007-2018\nMass Balance (m w.e. a$^{-1}$)',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0.8,0.88)
plt.grid()
fig.legend(fontsize=14,bbox_to_anchor=(1.5,0.94))
plt.tight_layout()
plt.savefig(os.path.join(output_stageII_params,'mb_vs_snowlinescore_scatterplot.pdf'),bbox_inches='tight')


# =============================================================================
# # Compare final parameters selected for the reference model to those from
# the alternative debris and accumulation models (fig 5.15 from KR thesis)
# =============================================================================
if compare_alt_model_params == True:
    # debris-free model
    NO_DEBRIS_params = np.loadtxt(debris_free_finalparams,skiprows=1,delimiter=',') 
    aice_final_NO_DEBRIS = NO_DEBRIS_params[:,1]
    asnow_final_NO_DEBRIS = NO_DEBRIS_params[:,2]
    MF_final_NO_DEBRIS = NO_DEBRIS_params[:,3]

    # Rounce et al. (2021) debris model
    rouncedebris_params = np.loadtxt(rounce_debris_finalparams,skiprows=1,delimiter=',') 
    aice_final_rouncedebris = rouncedebris_params[:,1]
    asnow_final_rouncedebris = rouncedebris_params[:,2]
    MF_final_rouncedebris = rouncedebris_params[:,3]

    # model with uncorrected accumulation:
    UNCORRECTED_ACC_params = np.loadtxt(uncorrected_acc_finalparams,skiprows=1,delimiter=',') 
    aice_final_uncorrected_acc = UNCORRECTED_ACC_params[:,1]
    asnow_final_uncorrected_acc = UNCORRECTED_ACC_params[:,2]
    MF_final_uncorrected_acc = UNCORRECTED_ACC_params[:,3]

    # model with AWS bias correction
    AWS_BIAS_params = np.loadtxt(aws_bc_finalparams,skiprows=1,delimiter=',') 
    aice_final_aws_bias = AWS_BIAS_params[:,1]
    asnow_final_aws_bias = AWS_BIAS_params[:,2]
    MF_final_aws_bias = AWS_BIAS_params[:,3]

    fig, axs = plt.subplots(2, 3, figsize=(7,3.8))

    axs[0,0].text(0.2e-6,26,'a)',fontsize=10,weight='bold')
    axs[0,1].text(0.1e-6,26,'b)',fontsize=10,weight='bold')
    axs[0,2].text(0.2e-4,26,'c)',fontsize=10,weight='bold')
    axs[1,0].text(0.2e-6,26,'d)',fontsize=10,weight='bold')
    axs[1,1].text(0.1e-6,26,'e)',fontsize=10,weight='bold')
    axs[1,2].text(0.2e-4,26,'f)',fontsize=10,weight='bold')

    axs[0,0].grid(zorder=-20)
    axs[0,0].hist(aice_final,color='grey',alpha=0.5,zorder=10,label='Reference model')
    axs[0,0].hist(aice_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[0,0].hist(aice_final_NO_DEBRIS,color='orange',alpha=0.5,zorder=10,label='Debris-free model')
    axs[0,0].hist(aice_final_NO_DEBRIS,color='orange',histtype='step',linewidth=3,zorder=10)
    axs[0,0].hist(aice_final_rouncedebris,color='maroon',alpha=0.5,zorder=10,label='Rounce et al. (2021) debris model')
    axs[0,0].hist(aice_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
    axs[0,0].set_ylabel('Frequency',fontsize=10)
    axs[0,0].set_xlabel('a$_{ice}$',fontsize=11)
    axs[0,0].set_xlim(0,np.max(aice_final)+0.4e-6)
    axs[0,0].set_xticks(np.arange(0,4.01e-6,1e-6))
    axs[0,0].tick_params(axis='both',labelsize=10)
    axs[0,0].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[0,0].set_ylim(0,30)
    axs[0,0].set_yticks(np.arange(0,31,5))

    axs[0,1].grid(zorder=-20)
    axs[0,1].hist(asnow_final,color='grey',alpha=0.5,zorder=10)
    axs[0,1].hist(asnow_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[0,1].hist(asnow_final_NO_DEBRIS,color='orange',alpha=0.5,zorder=10)
    axs[0,1].hist(asnow_final_NO_DEBRIS,color='orange',histtype='step',linewidth=3,zorder=10)
    axs[0,1].hist(asnow_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
    axs[0,1].hist(asnow_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
    #axs[0,1].set_ylabel('Frequency',fontsize=10)
    axs[0,1].set_xlabel('a$_{snow}$',fontsize=11)
    axs[0,1].set_xlim(0,np.max(asnow_final)+0.3e-6)
    axs[0,1].set_xticks(np.arange(0,2.01e-6,0.5e-6))
    axs[0,1].tick_params(axis='both',labelsize=10)
    axs[0,1].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[0,1].set_ylim(0,30)
    axs[0,1].set_yticks(np.arange(0,31,5))

    axs[0,2].grid(zorder=-20)
    axs[0,2].hist(MF_final,color='grey',alpha=0.5,zorder=10)
    axs[0,2].hist(MF_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[0,2].hist(MF_final_NO_DEBRIS,color='orange',alpha=0.5,zorder=10)
    axs[0,2].hist(MF_final_NO_DEBRIS,color='orange',histtype='step',linewidth=3,zorder=10)
    axs[0,2].hist(MF_final_rouncedebris,color='maroon',alpha=0.5,zorder=10)
    axs[0,2].hist(MF_final_rouncedebris,color='maroon',histtype='step',linewidth=3,zorder=10)
    #axs[0,2].set_ylabel('Frequency',fontsize=10)
    axs[0,2].set_xlabel('$MF$',fontsize=11)
    axs[0,2].set_xlim(0,np.max(MF_final)+0.3e-4)
    axs[0,2].set_xticks(np.arange(0,5.01e-4,1e-4))
    axs[0,2].tick_params(axis='both',labelsize=10)
    axs[0,2].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[0,2].set_ylim(0,30)
    axs[0,2].set_yticks(np.arange(0,31,5))

    axs[1,0].grid(zorder=-20)
    axs[1,0].hist(aice_final,color='grey',alpha=0.5,zorder=10)
    axs[1,0].hist(aice_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[1,0].hist(aice_final_uncorrected_acc,color='cornflowerblue',alpha=0.5,zorder=10,label='Model with uncorrected accumulation')
    axs[1,0].hist(aice_final_uncorrected_acc,color='cornflowerblue',histtype='step',linewidth=3,zorder=10)
    axs[1,0].hist(aice_final_aws_bias,color='indigo',alpha=0.5,zorder=10,label='Model with AWS bias correction')
    axs[1,0].hist(aice_final_aws_bias,color='indigo',histtype='step',linewidth=3,zorder=10)
    axs[1,0].set_ylabel('Frequency',fontsize=10)
    axs[1,0].set_xlabel('a$_{ice}$',fontsize=11)
    axs[1,0].set_xlim(0,np.max(aice_final)+0.4e-6)
    axs[1,0].set_xticks(np.arange(0,4.01e-6,1e-6))
    axs[1,0].tick_params(axis='both',labelsize=10)
    axs[1,0].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[1,0].set_ylim(0,30)
    axs[1,0].set_yticks(np.arange(0,31,5))

    axs[1,1].grid(zorder=-20)
    axs[1,1].hist(asnow_final,color='grey',alpha=0.5,zorder=10)
    axs[1,1].hist(asnow_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[1,1].hist(asnow_final_uncorrected_acc,color='cornflowerblue',alpha=0.5,zorder=10)
    axs[1,1].hist(asnow_final_uncorrected_acc,color='cornflowerblue',histtype='step',linewidth=3,zorder=10)
    axs[1,1].hist(asnow_final_aws_bias,color='indigo',alpha=0.5,zorder=10)
    axs[1,1].hist(asnow_final_aws_bias,color='indigo',histtype='step',linewidth=3,zorder=10)
    #axs[1,1].set_ylabel('Frequency',fontsize=10)
    axs[1,1].set_xlabel('a$_{snow}$',fontsize=11)
    axs[1,1].set_xlim(0,np.max(asnow_final)+0.3e-6)
    axs[1,1].set_xticks(np.arange(0,2.01e-6,0.5e-6))
    axs[1,1].tick_params(axis='both',labelsize=10)
    axs[1,1].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[1,1].set_ylim(0,30)
    axs[1,1].set_yticks(np.arange(0,31,5))

    axs[1,2].grid(zorder=-20)
    axs[1,2].hist(MF_final,color='grey',alpha=0.5,zorder=10)
    axs[1,2].hist(MF_final,color='grey',histtype='step',linewidth=3,zorder=10)
    axs[1,2].hist(MF_final_uncorrected_acc,color='cornflowerblue',alpha=0.5,zorder=10)
    axs[1,2].hist(MF_final_uncorrected_acc,color='cornflowerblue',histtype='step',linewidth=3,zorder=10)
    axs[1,2].hist(MF_final_aws_bias,color='indigo',alpha=0.5,zorder=10)
    axs[1,2].hist(MF_final_aws_bias,color='indigo',histtype='step',linewidth=3,zorder=10)
    #axs[1,2].set_ylabel('Frequency',fontsize=10)
    axs[1,2].set_xlabel('$MF$',fontsize=11)
    axs[1,2].set_xlim(0,np.max(MF_final)+0.3e-4)
    axs[1,2].set_xticks(np.arange(0,5.01e-4,1e-4))
    axs[1,2].tick_params(axis='both',labelsize=10)
    axs[1,2].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    axs[1,2].set_ylim(0,30)
    axs[1,2].set_yticks(np.arange(0,31,5))


    fig.legend(fontsize=10,ncol=2,bbox_to_anchor=(0.97,0.94),loc='lower right')

    fig.tight_layout()
    fig.savefig(os.path.join(output_stageII_params,'all_models_tuned_params.pdf'),bbox_inches='tight')

else:
    pass







