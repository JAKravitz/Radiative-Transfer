#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot EAP optics data

@author: jkravz311
"""
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

with open('/Users/jkravz311/Desktop/EAP_optics.p', 'rb') as fp:
    data = pickle.load(fp)

# lambda
l = np.arange(400,905,5)

#%%

fig, axs = plt.subplots(4,4,figsize=(20,20))
axs = axs.ravel()
count = 0
for c in data:
    for phyto in data[c]:
        for sname in data[c][phyto]:
            s = data[c][phyto][sname]['bb'] 
            for deff in s.index:
                axs[count].plot(l,s.loc[deff,:])
    axs[count].set_title(c)
    count = count + 1
            
#%%
import statsmodels.api as sm


colors = [(0, 1, 0), (0, .7, 0), (0, .3, 0)]
cmap_name = 'my_list'
n_bin=3


fig, axs = plt.subplots(4,4,figsize=(20,20))
axs = axs.ravel()
count = 0
for c in data:
    for phyto in data[c]:
        for sname in data[c][phyto]:
            s = data[c][phyto][sname]['bb'] 
            for deff in s.index:
                axs[count].plot(l,s.loc[deff,:])
    axs[count].set_title(c)
    count = count + 1

#%%
cyano = data['Cyanophyceae']['S. elongatus']['0.12_6.50_1.09']['a']


