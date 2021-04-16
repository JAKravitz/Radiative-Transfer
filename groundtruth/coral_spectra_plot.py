#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 21:44:38 2021

@author: jkravz311
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import statsmodels.api as sm

path = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/corals/reef_spectra.csv'
data = pd.read_csv(path)

coral = data[data['Comment'] == 'Stony coral']
coralspec = coral.iloc[:,4:-100].astype('float')

def kmeans_clust (data,k):
    # k means clustering
    kmeans = KMeans(n_clusters=k)
    kmeans = kmeans.fit(data)
    labels = kmeans.predict(data)
    #centroids = kmeans.cluster_centers_
    return labels

labels = kmeans_clust(coralspec,8)

coralspec['labels'] = labels

grouped = coralspec.groupby('labels')
for k,group in grouped:
    group.iloc[:,:-1].T.plot(legend=False)
    

#%%
from __future__ import division

l = np.arange(400,801,1)
fig, axs = plt.subplots(5,7,figsize=(25,20))
axs = axs.ravel()
count = 0
grouped = data.groupby('Comment')
for k,group in grouped:
    print (k)
    if group.shape[0] > 10:
        res = sm.graphics.fboxplot(group.iloc[:,4:].values, l, wfactor=50, ax=axs[count])
        axs[count].set_title(k)
        axs[count].set_ylim(0,3)
        count = count+1
    else:
        group.iloc[:,4:].T.plot(ax=axs[count],legend=False)
        #axs[count].plot(l, group.iloc[:,4:])
        axs[count].set_title(k)
        axs[count].set_ylim(0,3)
        count = count+1
fig.savefig('/Users/jkravz311/Desktop/coral_spectra.png',bbox_inches='tight',dpi=300)        
        
        
        
        