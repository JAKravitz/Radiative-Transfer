#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot EAP optics data

@author: jkravz311
"""
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

with open('/Users/jkravz311/Desktop/EAP_optics.json', 'r') as f:
    loaddata = json.load(f)
    data = json.loads(loaddata)

tert = data['Chlorophyceae']['D. tertiolecta1']
theta = data['Chlorophyceae']['G. theta']
data['Cryptophyceae']['G. theta'] = theta
data['Chlorophyceae'] = {'D. tertiolecta1':tert} 


#%%

# lambda
l = np.arange(400,905,5)

# absorption
fig, axs = plt.subplots(4,4,figsize=(20,20))
axs = axs.ravel()
count = 0
for c in data:
    for phyto in data[c]:
        for sname in data[c][phyto]:
            s = data[c][phyto][sname]['bbtilde'] 
            for deff in s:
                axs[count].plot(l,deff)
    count = count + 1
            
        