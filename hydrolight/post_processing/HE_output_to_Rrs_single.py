#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick look at hydrolight radiance output

@author: jkravz311
"""
#%%
import os
import re
import pandas as pd
import numpy as np

fname = '/Users/jkravz311/git_projects/Radiative-Transfer/hydrolight/test/Msunglint1.txt'

# flist = os.listdir(path)


wl = []
rrs = []
rrsfile = pd.DataFrame()
ed = []
lu = []
lw = []
with open(fname, 'r+') as f:
    for i,line in enumerate(f):

        # if i == 0:
        #     info = line
        #     title = re.findall(r'I_.*txt', info)

        if 434 < i < 535:
            wl.append(line.split(' ')[1])
            rrs.append(line.split(' ')[4])
            lw.append(line.split(' ')[10])
            lu.append(line.split(' ')[13].strip())


lu = np.array(lu,dtype=float)
lw = np.array(lw,dtype=float)
glint = lu - lw
rrsfile['lu'] = lu
rrsfile['lw'] = lw
rrsfile['glint'] = glint
rrsfile.index = np.array(wl,dtype=float)

rrsfile.plot()

