# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:45:10 2023

Generates N random combinations of a_snow, a_ice, MF, and saves them to a csv to be used in
model tuning.

Parameters are randomly selected from independent normal distributions truncated at zero based on the 
mean and std. dev. of values of a_snow, a_ice, MF found in literature (see Young et al. 2021)

@author: katierobinson
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

np.random.seed(6) # Set seed to make the parameter set reproducible

# N is the number of parameters to be generated.
N = 10000

# mean and std. dev. of values found in literature (see Young et al. 2021)
aice_mean = 0.000003396
aice_std = 0.00000438

asnow_mean = 0.000001546
asnow_std = 0.00000085

MF_mean = 0.0002707
MF_std = 0.0001632

aice_list = []
asnow_list = []
MF_list = []
for i in range(0,N):

    aice = np.random.normal(aice_mean, aice_std) 
    while aice < 0:
        aice = np.random.normal(aice_mean, aice_std)
    aice_list.append(aice)
    
    asnow = np.random.normal(asnow_mean,asnow_std)
    while asnow < 0:
        asnow = np.random.normal(asnow_mean,asnow_std)
    asnow_list.append(asnow)
        
    MF = np.random.normal(MF_mean,MF_std)
    while MF < 0:
        MF = np.random.normal(MF_mean,MF_std)
    MF_list.append(MF)

d = {'aice': aice_list, 'asnow': asnow_list, 'MF': MF_list}
df = pd.DataFrame(data=d)

df.to_csv('Tuning_Params.csv')

# =============================================================================
# Plot param distributions
# =============================================================================

aice_bins = np.linspace(0,2.1e-5,int(np.sqrt(N)))
asnow_bins = np.linspace(0,4.7e-6,int(np.sqrt(N)))
MF_bins = np.linspace(0,8.9e-4,int(np.sqrt(N)))

fig, (ax) = plt.subplots(ncols=3,figsize=(7,2.5),gridspec_kw={"width_ratios":[1,1,1]},sharey=True)
fig.subplots_adjust(wspace=0.3)

plt.subplot(1,3,1)
plt.title('a$_{ice}$',fontsize=12)
plt.hist(aice_list,aice_bins,rwidth=1,color='lightgrey',label='All simulations')
plt.ylim(0,300)
plt.ylabel('Frequency',fontsize=12)
plt.vlines(x=0.000003396,ymin=0,ymax=2500,linestyle='--')
x = np.linspace(np.min(aice_list), np.max(aice_list), 100)
p = norm.pdf(x, 0.000003396, 0.00000438)
plt.plot(x, p/np.max(p)*np.max(np.histogram(aice_list,aice_bins)[0]), 'k', linewidth=2,label='Normal distribution')
plt.margins(x=0)
plt.xticks(np.arange(0,2.1e-5,0.5e-5))
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

plt.subplot(1,3,2)
plt.title('a$_{snow}$',fontsize=12)
plt.hist(asnow_list,asnow_bins,rwidth=1,color='lightgrey')
plt.ylim(0,300)
plt.vlines(x=0.000001546,ymin=0,ymax=2500,linestyle='--')
x = np.linspace(np.min(asnow_list), np.max(asnow_list), 100)
p = norm.pdf(x, 0.000001546,0.00000085)
plt.plot(x, p/np.max(p)*np.max(np.histogram(asnow_list,asnow_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)
plt.xticks(np.arange(0,4.5e-6,1e-6))
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

plt.subplot(1,3,3)
plt.title('$MF$',fontsize=12)
plt.hist(MF_list,MF_bins,rwidth=1,color='lightgrey')
plt.ylim(0,300)
plt.vlines(x=0.0002707,ymin=0,ymax=2500,linestyle='--')
x = np.linspace(np.min(MF_list), np.max(MF_list), 100)
p = norm.pdf(x, 0.0002707,0.0001632)
plt.plot(x, p/np.max(p)*np.max(np.histogram(MF_list,MF_bins)[0]), 'k', linewidth=2)
plt.margins(x=0)
plt.xticks(np.arange(0,9.1e-4,2e-4))
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

plt.tight_layout()
#plt.savefig('F:\Mass Balance Model\Kaskawulsh-Mass-Balance\Tuning\MeltModelParams_histograms.pdf',bbox_inches='tight')

