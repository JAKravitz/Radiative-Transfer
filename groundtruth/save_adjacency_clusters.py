#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Working script...
Organizing aster spectral library  for adjacency

kravitz, 2021
"""
#%%
import pandas as pd
import numpy as np

datapath = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/aster_spectral_library/jpl_aster_data.csv'
data = pd.read_csv(datapath)
meta = data.filter(regex='\D')
spectra = data.filter(regex='\d')

#%% water

snow = data[data.Type == 'Water']
snow = snow.drop([64,157])


#%% manmade

man = data[data.Type == 'manmade']

# concrete spectra
concrete = man[man.Class == 'Concrete']

# construction material spectra
construction = man[man.Class == 'General Construction Material']
construction = construction.loc[[175,225,364,467,562],:]

# Road spectra
road = man[man.Class == 'Road']
roadspec = road.filter(regex='\d')
# first 20nm are bunk, make 0
bunkcols = roadspec.iloc[:,:21]
bunk = np.where(bunkcols > 0, np.nan, np.nan)
bunk = pd.DataFrame(bunk,index=bunkcols.index,columns=bunkcols.columns)
roadspec = roadspec.iloc[:,21:]
roadmeta = road.filter(regex ='\D')
road = pd.concat([roadmeta,bunk,roadspec],axis=1)

# roofing material spectra
roof = man[man.Class == 'Roofing Material']

#%% non photosynthetic vegetation

nonveg = data[data.Type=='non photosynthetic vegetation']
nonveg = nonveg[(nonveg.Class != 'flowers') & (nonveg.Class != 'leaves') & (nonveg.Class != 'lichen')]

#%% vegetation

data = data.replace(to_replace='Vegetation',value='vegetation')
veg = data[data.Type == 'vegetation']

# shrub
shrub = veg[veg.Class == 'Shrub']

# tree
tree = veg[veg.Class == 'Tree']

# grass
grass = veg[veg.Class == 'grass']

#%% soil

soil = data[data.Type == 'soil']
soil = soil[(soil.Class != 'Utisol') & (soil.Class != 'Vertisol')]

#%% check linear mixing model

import matplotlib.pyplot as plt

soil1 = soil.sample().filter(regex='\d').values * .5
tree1 = tree.sample().filter(regex='\d').values * .3
road1 = road.sample().filter(regex='\d').values * .2
mix = np.nansum(np.stack((soil1[0],tree1[0],road1[0])), axis=0)
#mix = soil1[0]+tree1[0]+road1[0]

wl = np.arange(400,900,1)

fig,ax=plt.subplots()
ax.plot(wl,soil1[0]*2,label='orig soil')
ax.plot(wl,soil1[0],label='soil')
ax.plot(wl,tree1[0],label='tree')
ax.plot(wl,road1[0],label='road')
ax.plot(wl,mix,label='mix')
ax.legend()

#%% Save
snow.to_csv('/Users/jkravz311/Desktop/snow_adjacency.csv',index=False)
concrete.to_csv('/Users/jkravz311/Desktop/concrete_adjacency.csv',index=False)
construction.to_csv('/Users/jkravz311/Desktop/construction_adjacency.csv',index=False)
road.to_csv('/Users/jkravz311/Desktop/road_adjacency.csv',index=False)
roof.to_csv('/Users/jkravz311/Desktop/roof_adjacency.csv',index=False)
nonveg.to_csv('/Users/jkravz311/Desktop/non_photo_veg_adjacency.csv',index=False)
shrub.to_csv('/Users/jkravz311/Desktop/shrub_adjacency.csv',index=False)
tree.to_csv('/Users/jkravz311/Desktop/tree_adjacency.csv',index=False)
grass.to_csv('/Users/jkravz311/Desktop/grass_adjacency.csv',index=False)
soil.to_csv('/Users/jkravz311/Desktop/soil_adjacency.csv',index=False)

#%%
# clusters = soil.groupby('Class')
# clustDict = {}
# for key, group in clusters:
#     cluster = group.filter(regex='\d')
#     clustDict[key] = cluster
#     cluster.T.plot(legend=False)
#     print (key)



    