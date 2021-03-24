#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Working script

@author: jkravz311
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('/Users/jkravz311/Desktop/phyto_imag_refr_idx.csv')
#spectra = data.filter(regex='^[0-9]')
l = np.arange(400,801,1)
grouped = data.groupby('Class')

fig, axs = plt.subplots(4,3,figsize=(20,20))
axs = axs.ravel()
count = 0
for k, group in grouped:
    spectra = group.filter(regex='^[0-9]')
    spectra.T.plot(ax=axs[count])
    axs[count].set_title(k)
    axs[count].set_ylim(0,.15)
    count = count +1

fig.savefig('/Users/jkravz311/Desktop/phyto_aphy.png',bbox_inches='tight',dpi=300)

