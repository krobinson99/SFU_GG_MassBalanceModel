# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 11:33:55 2023

Determines which simulations pass stage I of tuning. To pass stage 1, simulations must
1) have a value of aice >= asnow
2) have a modelled 2007-2018 net mass balance within 3 std. dev of the geodetic mass balance (-0.46 m w.e. a-1)

Inputs: 
- params_file: a csv file containing the values of aice, asnow, and MF used for tuning
- modelled_mb: a csv file containing the modelled mass balance from each of the parameter combinations in the params_file
- optional: csv files with modelled mass balance results for the debris and accumulation sensitivity tests

Outputs:
- a csv file containing the values of aice, asnow, and MF from sims that pass stage 1
- a csv file containing the sim ID and modelled mass balance from sims that pass stage 1

This script also plots figures 5.1, 5.2, 5.6 and 5.10 from KR's thesis

@author: katierobinson
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import os

# =============================================================================
# Paths to melt model parameters and modelled 2007-2018 net mass balances
# =============================================================================
model_name = 'REF_MODEL' # string assigned to this run of the model to help keep track of results
params_file = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/meltmodel_params/tuning_params_0-9999.csv' # Path to csv file containing aice, asnow, MF used for tuning
modelled_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/refmodel_mb_sim_0-9999.csv' #  csv file containing the net mb for each sim (calculated using calulate_mb.py)
output_stage1_params = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results' # directory where params that pass stage 1 should be stored

# optional: compare with modelled mb's from the alternative debris and accumulation tests (Robinson et al. 2024)
compare_alt_debris_models = True
rounce_debris_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/nodebris_mb_sim_0-9999.csv'
debris_free_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/rouncedebris_mb_sim_0-9999.csv'

compare_alt_accumulation_models = True
uncorrected_acc_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/uncorracc_mb_sim_0-9999.csv'
aws_bc_mb = '/home/krobin/projects/def-gflowers/krobin/MBM/Kaskawulsh/tuning_results/awsbias_mb_sim_0-9999.csv'

# =============================================================================
# Load melt model parameters and modelled mass balances
# =============================================================================
params = np.loadtxt(params_file,skiprows=1,delimiter=',') 
aice = params[:,1]
asnow = params[:,2]
MF = params[:,3]

mb = np.loadtxt(modelled_mb,delimiter=',') 
simID = np.arange(0,len(mb))  

# =============================================================================
# Remove sims where aice < asnow
# =============================================================================
mb_aice_geq_asnow = np.delete(mb,np.where(aice<asnow))
simID_aice_geq_asnow = np.delete(simID,np.where(aice<asnow))

# =============================================================================
# Define the mass balance target
# =============================================================================
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

# Save params that pass tuning stage I (aice >= asnow and meets MB target) in new csv:

d = {'aice': list(aice_passing), 'asnow': list(asnow_passing), 'MF': list(MF_passing)}
df = pd.DataFrame(data=d)
df.to_csv(os.path.join(output_stage1_params,'tuning_stageI_params_' + model_name + '.csv'))

d2 = {'simID': simID_passing, 'MB': mb_passing}
df2 = pd.DataFrame(data=d2)
df2.to_csv(os.path.join(output_stage1_params,'tuning_stageI_mb_' + model_name + '.csv'))

# =============================================================================
# Plot histogram of params that were tested (fig 5.1 from KR thesis)
# =============================================================================

# Get simIDs and param values for all sims that pass MB target (regardless of aice < or > asnow)
simID_passing_allsims = simID[np.where((mb >= (target_mb-(3*tagret_mb_stddev))) & (mb <= (target_mb+(3*tagret_mb_stddev))))]
mb_passing_allsims = mb[np.where((mb >= (target_mb-(3*tagret_mb_stddev))) & (mb <= (target_mb+(3*tagret_mb_stddev))))]
aice_passing_allsims = aice[simID_passing_allsims]
asnow_passing_allsims = asnow[simID_passing_allsims]
MF_passing_allsims = MF[simID_passing_allsims]

N = int(np.sqrt(len(mb)))

aice_bins = np.linspace(0,2.1e-5,N)
asnow_bins = np.linspace(0,4.7e-6,N)
MF_bins = np.linspace(0,8.9e-4,N)

aicestd = 0.00000438
asnowstd = 0.00000085
MFstd = 0.0001632

fig, (ax) = plt.subplots(ncols=3,figsize=(9,3.5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=14)
plt.hist(aice,aice_bins,rwidth=1,color='lightgrey',label='All simulations')
plt.hist(aice_passing_allsims,aice_bins,rwidth=1,color='red',zorder=5,alpha=0.5,label='$\dot{B}_{mod}$ = $\dot{B}_{obs}$ $\pm$ 3$\sigma$')
plt.hist(aice_passing,aice_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5,label='$\dot{B}_{mod}$ = $\dot{B}_{obs}$ $\pm$ 3$\sigma$\nand a$_{ice}$ $\geq$ a$_{snow}$')
plt.hist(aice_passing_allsims,aice_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(aice_passing,aice_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)

plt.ylim(0,300)
plt.ylabel('Frequency',fontsize=14)
plt.vlines(x=0.000003396,ymin=0,ymax=2500,linestyle='--')
#plt.xticks([0,0.000003396,0.000003396+(1*aicestd),0.000003396+(2*aicestd),0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','7.8x10$^{-6}$','1.2x10$^{-5}$','1.6x10$^{-5}$'],rotation=45)
plt.xticks([0,0.000003396,0.000003396+(3*aicestd)],['0','3.4x10$^{-6}$','1.6x10$^{-5}$'],rotation=0)
x = np.linspace(np.min(aice), np.max(aice), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=14)
plt.hist(asnow,asnow_bins,rwidth=1,color='lightgrey')
plt.hist(asnow_passing_allsims,asnow_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(asnow_passing,asnow_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
plt.hist(asnow_passing_allsims,asnow_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(asnow_passing,asnow_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,300)
plt.vlines(x=0.000001546,ymin=0,ymax=2500,linestyle='--')
plt.xticks([0,0.000001546,0.000001546+(3*asnowstd)],['0','1.5x10$^{-6}$','4.1x10$^{-6}$'])
x = np.linspace(np.min(asnow), np.max(asnow), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=14)
plt.hist(MF,MF_bins,rwidth=1,color='lightgrey')
plt.hist(MF_passing_allsims,MF_bins,rwidth=1,color='red',zorder=5,alpha=0.5)
plt.hist(MF_passing,MF_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5)
plt.hist(MF_passing_allsims,MF_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=1.5,zorder=10)
plt.hist(MF_passing,MF_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=1.5,zorder=10)
plt.ylim(0,300)
plt.vlines(x=0.0002707,ymin=0,ymax=2500,linestyle='--')
plt.xticks([0,0.0002707,0.0002707+(3*MFstd)],['0','2.7x10$^{-4}$','7.6x10$^{-4}$'])
x = np.linspace(np.min(MF), np.max(MF), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)

#handles, labels = ax.get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
fig.legend(fontsize=14,bbox_to_anchor=(1.18,0.8), ncol=1, borderaxespad=0.19)
plt.tight_layout()
plt.savefig(os.path.join(output_stage1_params,'meltmodelparams_priordistributions.pdf'),bbox_inches='tight')

# =============================================================================
# Plot distribution modelled mass balances from stage 1 (fig 5.2 from KR thesis)
# =============================================================================
mb_aice_less_asnow = np.delete(mb,np.where(aice>=asnow))
simID_aice_less_asnow = np.delete(simID,np.where(aice>=asnow))

target_mb = -0.46
tagret_mb_stddev = 0.17
mb_discarded = mb_aice_less_asnow[np.where((mb_aice_less_asnow >= (target_mb-(3*tagret_mb_stddev))) & (mb_aice_less_asnow <= (target_mb+(3*tagret_mb_stddev))))]
print('the number of discarded simulations with the correct MB is ',len(mb_discarded))
mb_bins =  np.linspace(min(mb),max(mb),100)

plt.figure(figsize=(8.4,3.5))
plt.text(-11.5,310,'Site-specific reference model',weight='bold',fontsize=14)
#plt.title('All runs\nTotal = ' + str(len(mb)),fontsize=14)
plt.hist(mb,mb_bins,rwidth=1,color='red',zorder=5,alpha=0.5,label='All simulations\nTotal = 10000')
plt.hist(mb_aice_geq_asnow,mb_bins,rwidth=1,color='mediumblue',zorder=5,alpha=0.5,label='Simulations where a$_{ice}$ $\geq$ a$_{snow}$\nTotal = ' + str(len(mb_aice_geq_asnow)))
plt.hist(mb,mb_bins, histtype='step', stacked=True, fill=False,color='red',linewidth=3,zorder=10)
plt.hist(mb_aice_geq_asnow,mb_bins, histtype='step', stacked=True, fill=False,color='mediumblue',linewidth=3,zorder=10)
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.xticks(np.arange(-11.5,1.1,0.5),np.arange(-11.5,1.1,0.5),fontsize=14,rotation=45)
plt.yticks(fontsize=14)
plt.ylim(0,300)
plt.xlim(-12,1.5)
plt.grid()
plt.legend(fontsize=14)
plt.vlines(-0.46-(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
plt.vlines(-0.46+(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
plt.hlines(1,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
plt.hlines(299,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
plt.tight_layout()
plt.savefig(os.path.join(output_stage1_params,'modelled_massbalance_distribution.pdf'),bbox_inches='tight')

# =============================================================================
# Plot mass balance results that fall under the MB normal distribution (scaled to include 100 sims)
# =============================================================================
mb_normaldist_bins =  np.linspace(target_mb-(3*tagret_mb_stddev),target_mb+(3*tagret_mb_stddev),int(np.sqrt(len(mb_passing))))
mb_normaldist_bins =  np.linspace(target_mb-(3*tagret_mb_stddev),target_mb+(3*tagret_mb_stddev),24)
bin_centers = np.arange(target_mb-(3*tagret_mb_stddev)+(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,target_mb+(3.1*tagret_mb_stddev)-(mb_normaldist_bins[1] - mb_normaldist_bins[0])/2,(mb_normaldist_bins[1] - mb_normaldist_bins[0]))

plt.figure(figsize=(12,5))
plt.grid(zorder=20)
plt.title('Target MB $\pm$ 3$\sigma$\nTotal = ' + str(len(mb_passing)),fontsize=14)
plt.hist(mb_passing,mb_normaldist_bins,rwidth=0.9,color='mediumblue',zorder=10)
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
plt.yticks(fontsize=14)
#plt.xticks(np.arange(-0.98,0.12,0.04),np.round(np.arange(-0.98,0.12,0.04),2),fontsize=14,rotation=45)
plt.xticks(bin_centers,np.round(bin_centers,2),fontsize=14,rotation=45)
x = np.arange(-0.98,0.11,0.02) # aligns with the center of each bin
p = norm.pdf(x, target_mb,tagret_mb_stddev)
#plt.ylim(0,30)
plt.plot(x, p/np.max(p)*np.histogram(mb_aice_geq_asnow,mb_normaldist_bins)[0][11], 'k', linewidth=4,zorder=30)
plt.tight_layout()

# x[np.where((p/np.max(p)*11)>=1)] # 11 is the number of sims in the bin centered on target_mb
# approx range of mb_vals inside normal distribution is -1 to 0.1

if compare_alt_debris_models == True:
    # =============================================================================
    # Compare mass balance distributions from ref model with the various debris treatments
    # =============================================================================

    # Load debris free results:
    mb_nodebris = np.loadtxt(debris_free_mb,delimiter=',')
    mb_aice_geq_asnow_nodebris = np.delete(mb_nodebris,np.where(aice<asnow))

    # Load rounce debris results:
    mb_rouncedeb = np.loadtxt(rounce_debris_mb,delimiter=',')
    mb_aice_geq_asnow_rouncedeb = np.delete(mb_rouncedeb,np.where(aice<asnow))

    # =============================================================================
    # Mass balance distribution of ref model vs alternative debris models (fig. 5.6 from KR thesis)
    # =============================================================================

    N = 100
    mb_bins =  np.linspace(min(mb),max(mb),N)

    plt.figure(figsize=(8.4,7))
    plt.subplot(2,1,1)
    #plt.title('All runs\nTotal = ' + str(len(mb)),fontsize=14)
    plt.hist(mb,mb_bins,rwidth=1,color='grey',zorder=5,alpha=0.5,label='Reference Model')
    plt.hist(mb_nodebris,mb_bins,rwidth=1,color='orange',zorder=5,alpha=0.5,label='Debris-free model')
    plt.hist(mb_rouncedeb,mb_bins,rwidth=1,color='maroon',zorder=5,alpha=0.5,label='Rounce et al. (2021)\ndebris model')
    plt.hist(mb,mb_bins, histtype='step', stacked=True, fill=False,color='grey',linewidth=3,zorder=10)
    plt.hist(mb_nodebris,mb_bins, histtype='step', stacked=True, fill=False,color='orange',linewidth=3,zorder=10)
    plt.hist(mb_rouncedeb,mb_bins, histtype='step', stacked=True, fill=False,color='maroon',linewidth=3,zorder=10)
    plt.ylabel('Frequency',fontsize=14)
    plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
    plt.xticks(np.arange(-11.5,1.1,0.5),np.arange(-11.5,1.1,0.5),fontsize=14,rotation=45)
    plt.yticks(fontsize=14)
    plt.ylim(0,300)
    plt.xlim(-12,1.5)
    plt.grid()
    plt.legend(fontsize=14)
    plt.vlines(-0.46-(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.vlines(-0.46+(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(1,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(299,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.text(-11.5,310,'a) All simulations (total = 10,000)',fontsize=14,weight='bold')

    N = int(np.sqrt(8394))
    mb_bins =  np.linspace(min(mb),max(mb),N)

    plt.subplot(2,1,2)
    #plt.title('All runs\nTotal = ' + str(len(mb)),fontsize=14)
    plt.hist(mb_aice_geq_asnow,mb_bins,rwidth=1,color='grey',zorder=5,alpha=0.5,label='Reference Model')
    plt.hist(mb_aice_geq_asnow_nodebris,mb_bins,rwidth=1,color='orange',zorder=5,alpha=0.5,label='Debris-free model')
    plt.hist(mb_aice_geq_asnow_rouncedeb,mb_bins,rwidth=1,color='maroon',zorder=5,alpha=0.5,label='Rounce et al. (2021)\ndebris model')
    plt.hist(mb_aice_geq_asnow,mb_bins, histtype='step', stacked=True, fill=False,color='grey',linewidth=3,zorder=10)
    plt.hist(mb_aice_geq_asnow_nodebris,mb_bins, histtype='step', stacked=True, fill=False,color='orange',linewidth=3,zorder=10)
    plt.hist(mb_aice_geq_asnow_rouncedeb,mb_bins, histtype='step', stacked=True, fill=False,color='maroon',linewidth=3,zorder=10)
    plt.ylabel('Frequency',fontsize=14)
    plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
    plt.xticks(np.arange(-11.5,1.1,0.5),np.arange(-11.5,1.1,0.5),fontsize=14,rotation=45)
    plt.yticks(fontsize=14)
    plt.ylim(0,300)
    plt.xlim(-12,1.5)
    plt.grid()
    plt.legend(fontsize=14)
    plt.vlines(-0.46-(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.vlines(-0.46+(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(1,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(299,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.text(-11.5,310,'b) Simulations where a$_{ice}$ $\geq$ a$_{snow}$ (total = 8394)',fontsize=14,weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_stage1_params,'debris_massbalance_distributions.pdf'),bbox_inches='tight')
else:
    pass


if compare_alt_accumulation_models == True:
    # =============================================================================
    # Compare mass balance distributions from ref model with the various accumulation treatments
    # =============================================================================

    # Load uncorrected accumulation results:
    mb_uncorracc = np.loadtxt(uncorrected_acc_mb,delimiter=',')
    mb_aice_geq_asnow_uncorracc = np.delete(mb_uncorracc,np.where(aice<asnow))

    # Load rounce debris results:
    mb_awsbias = np.loadtxt(aws_bc_mb,delimiter=',')
    mb_aice_geq_asnow_awsbias = np.delete(mb_awsbias,np.where(aice<asnow))

    # =============================================================================
    # Mass balance distribution of ref model vs alternative accumulation models (fig. 5.10 from KR thesis)
    # =============================================================================

    N = 100
    mb_bins =  np.linspace(min(mb),max(mb),N)

    plt.figure(figsize=(8.4,7))
    plt.subplot(2,1,1)
    #plt.title('All runs\nTotal = ' + str(len(mb)),fontsize=14)
    plt.hist(mb,mb_bins,rwidth=1,color='grey',zorder=5,alpha=0.5,label='Reference Model')
    plt.hist(mb_uncorracc,mb_bins,rwidth=1,color='cornflowerblue',zorder=5,alpha=0.5,label='Uncorrected\naccumulation')
    plt.hist(mb_awsbias,mb_bins,rwidth=1,color='indigo',zorder=5,alpha=0.5,label='ECCC AWS bias\ncorrection')
    plt.hist(mb,mb_bins, histtype='step', stacked=True, fill=False,color='grey',linewidth=3,zorder=10)
    plt.hist(mb_uncorracc,mb_bins, histtype='step', stacked=True, fill=False,color='cornflowerblue',linewidth=3,zorder=10)
    plt.hist(mb_awsbias,mb_bins, histtype='step', stacked=True, fill=False,color='indigo',linewidth=3,zorder=10)
    plt.ylabel('Frequency',fontsize=14)
    plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
    plt.xticks(np.arange(-11.5,1.1,0.5),np.arange(-11.5,1.1,0.5),fontsize=14,rotation=45)
    plt.yticks(fontsize=14)
    plt.ylim(0,300)
    plt.xlim(-12,1.5)
    plt.grid()
    plt.legend(fontsize=14)
    plt.vlines(-0.46-(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.vlines(-0.46+(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(1,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(299,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)

    plt.text(-11.5,310,'a) All simulations (total = 10,000)',fontsize=14,weight='bold')

    N = int(np.sqrt(8394))
    mb_bins =  np.linspace(min(mb),max(mb),N)

    plt.subplot(2,1,2)
    #plt.title('All runs\nTotal = ' + str(len(mb)),fontsize=14)
    plt.hist(mb_aice_geq_asnow,mb_bins,rwidth=1,color='grey',zorder=5,alpha=0.5,label='Reference Model')
    plt.hist(mb_aice_geq_asnow_uncorracc,mb_bins,rwidth=1,color='cornflowerblue',zorder=5,alpha=0.5,label='Uncorrected\naccumulation')
    plt.hist(mb_aice_geq_asnow_awsbias,mb_bins,rwidth=1,color='indigo',zorder=5,alpha=0.5,label='ECCC AWS bias\ncorrection')
    plt.hist(mb_aice_geq_asnow,mb_bins, histtype='step', stacked=True, fill=False,color='grey',linewidth=3,zorder=10)
    plt.hist(mb_aice_geq_asnow_uncorracc,mb_bins, histtype='step', stacked=True, fill=False,color='cornflowerblue',linewidth=3,zorder=10)
    plt.hist(mb_aice_geq_asnow_awsbias,mb_bins, histtype='step', stacked=True, fill=False,color='indigo',linewidth=3,zorder=10)
    plt.ylabel('Frequency',fontsize=14)
    plt.xlabel('2007-2018 Mass Balance (m w.e. $a^{-1}$)',fontsize=14)
    plt.xticks(np.arange(-11.5,1.1,0.5),np.arange(-11.5,1.1,0.5),fontsize=14,rotation=45)
    plt.yticks(fontsize=14)
    plt.ylim(0,300)
    plt.xlim(-12,1.5)
    plt.grid()
    plt.legend(fontsize=14)
    plt.vlines(-0.46-(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.vlines(-0.46+(3*0.17),ymin=0,ymax=300,linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(1,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.hlines(299,xmin=-0.46-(3*0.17),xmax=-0.46+(3*0.17),linewidth=3,linestyle='--',color='k',zorder=20)
    plt.text(-11.5,310,'b) Simulations where a$_{ice}$ $\geq$ a$_{snow}$ (total = 8394)',fontsize=14,weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_stage1_params,'accumulation_massbalance_distributions.pdf'),bbox_inches='tight')
else:
    pass