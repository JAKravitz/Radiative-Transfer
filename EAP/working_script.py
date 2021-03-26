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

data = pd.read_csv('/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/phyto_optics/in_vivo_phyto_abs.csv')
#spectra = data.filter(regex='^[0-9]')
l = np.arange(400,801,1)
grouped = data.groupby('Class')

fig, axs = plt.subplots(4,4,figsize=(20,20))
axs = axs.ravel()
count = 0
for k, group in grouped:
    print (k)
    spectra = group.filter(regex='^[0-9]')
    spectra.T.plot(ax=axs[count])
    axs[count].set_title(k)
    axs[count].set_ylim(0,.15)
    count = count +1

fig.savefig('/Users/jkravz311/Desktop/phyto_aphy.png',bbox_inches='tight',dpi=300)

#%%
from scipy.interpolate import griddata
import numpy as np
import pandas as pd

newdata = []

l = np.arange(.4, .905, .005) 
data = pd.read_csv('/Users/jkravz311/Desktop/W1.csv')
data = data.iloc[:,3:]
wv = data.columns.values.astype(float) / 1000

for k in data.index:
    s1 = data.iloc[k,:].values
    s2 = griddata(wv, s1, l, 'linear') # interpolate 
    s2 = np.where(s2 < 0, np.nan, s2) # replace negatives w/ nan
    s2 = np.where(s2 == 0, np.nan, s2) # replace zeros w/ nan
    s2[np.isnan(s2)] = min(s2)
    newdata.append(s2)

data2 = np.array(newdata)
data2 = pd.DataFrame(data2)
data2.columns = l

data2.to_csv('/Users/jkravz311/Desktop/W1_fin.csv')


#%%
s1 = data.iloc[0,:].values
s2 = griddata(wv, s1, l, 'linear') # interpolate 
s2 = np.where(s2 < 0, np.nan, s2) # replace negatives w/ nan
s2 = np.where(s2 == 0, np.nan, s2) # replace zeros w/ nan
s2[np.isnan(s2)] = min(s2)
newdata.append(s2)





# kshell_base = griddata(im_wv, im_a1, l, 'linear',)
# # replace negative/zeros with nans
# kshell_base = np.where(kshell_base < 0, np.nan, kshell_base)
# kshell_base = np.where(kshell_base == 0, np.nan, kshell_base)
# # replace nans with minimum value
# kshell_base[np.isnan(kshell_base)] = min(kshell_base)