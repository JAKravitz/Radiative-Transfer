#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Working script

@author: jkravz311
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from EAP_proc import EAP

#%%

optics = {'Bacillariophyceae': {},
          'Chlorophyceae': {},
          'Coscinodiscophyceae': {},
          'Cryptophyceae': {},
          'Cyanophyceae': {},
          'Dinophyceae': {},
          'Eustigmatophyceae': {},
          'Fragilariophyceae': {},
          'Pelagiophyceae': {},
          'Prasinophyceae': {},
          'Prymnesiophyceae': {},
          'Raphidophyceae': {},
          'Rhodophyceae': {}
          }

# phyto species
phytolist = ['D. tertiolecta1', 'D. tertiolecta2']
Class = 'Chlorophyceae'
# paths to imaginary and real refractive index
mf = '/Users/jkravz311/git_projects/Radiative-Transfer/EAP/501nm_extended_e1701000.mat'
astarpath = '/Users/jkravz311/git_projects/Radiative-Transfer/EAP/data/in_vivo_phyto_abs.csv'
batchinfo = pd.read_csv('/Users/jkravz311/git_projects/Radiative-Transfer/EAP/data/EAP_batch_V1.csv', index_col=0)
batchinfo = batchinfo.iloc[0:2,:]

## Optical parameters to vary ##

# Deff = effective diameter of phyto cell
# retrieved from batch info file
# varies by species
#
# Vs = relative chlorplast volume
# positive relationship with size
#
# Ci = intracellular chl density mg/m-3
# range changes by phyto type 
# maybe inverse relationship w/ size
# literature range: 0.1 - 12 mg/m-3
#
# nshell = real refractive index of outer core (shell)
# ranges 1.06-1.22

params = {'Vs1': np.arange(0.04,0.26,0.02),
          'Vs2': np.arange(0.2,0.4,0.02),
          'Vs3': np.arange(0.36,0.56,0.02),
          'Ci1': np.arange(.5e6,7e6,.5e6),
          'Ci2': np.arange(4e6,12e6,.5e6),
          'Ci3': np.arange(2e6,8e6,.5e6),
          'nshell1': np.arange(1.06,1.16,.01).astype('float16'),
          'nshell2': np.arange(1.11, 1.22,.01).astype('float16'),
          'nshell3': np.arange(1.1, 1.19,.01).astype('float16'),
          'ncore': np.arange(1.014, 1.04, 0.005)}

#%%
import json
import timeit as ti
from joblib import Parallel, delayed
import multiprocessing

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


num_cores = multiprocessing.cpu_count()

for i,phyto in enumerate(batchinfo.index):
    
    if i == 15:
        print ('$$$$$$$$ 25% DONE $$$$$$$$$')
    elif i == 30:
        print ('$$$$$$$$ 50% DONE $$$$$$$$$')
    elif i == 45:
        print ('$$$$$$$$ 75% DONE $$$$$$$$$')
    
    print ('####### {} #######'.format(phyto))
    
    optics[Class][phyto] = {}
    info = batchinfo.loc[phyto,:]
    VsF = np.random.choice(params[info.Vs], 3)
    CiF = np.random.choice(params[info.Ci], 3)
    nshellF = np.random.choice(params[info.nshell], 3)
    ncore = np.random.choice(params['ncore'], 1)
    Deff = np.arange(info.Dmin, info.Dmax, 1)
    
    for Vs in VsF:
        for ci in CiF:
            for n in nshellF:
                
                sname = '{:.2f}_{:.2f}_{:.2f}'.format(Vs, ci/1e6, n)
                
                # EAP run
                # standard
                tic = ti.default_timer()
                print ('------ {} ------'.format(sname))
                Qc, Sigma_c, c, Qb, Sigma_b, b, Qa, Sigma_a, a, Qbb, Sigma_bb, bb, bbtilde = EAP(phyto, mf, astarpath, Vs, ci, Deff, n, ncore)
                toc = ti.default_timer()
                runtime = toc - tic
                print (runtime)
                
                # in parallel
                # tic = ti.default_timer()
                # print ('------ {} ------'.format(sname))
                # Qc, Sigma_c, c, Qb, Sigma_b, b, Qa, Sigma_a, a, Qbb, Sigma_bb, bb, bbtilde = Parallel(n_jobs=num_cores)(delayed(EAP(phyto, mf, astarpath, Vs, ci, Deff, n, ncore)))
                # toc = ti.default_timer()
                # runtime = toc - tic
                # print (runtime)                
                
                # Store
                optics[Class][phyto][sname] = {'Qc':Qc,
                                               'Sigma_c':Sigma_c,
                                               'c':c,
                                               'Qb':Qb,
                                               'Sigma_b':Sigma_b,
                                               'b':b,
                                               'Qa':Qa,
                                               'Sigma_a':Sigma_a,
                                               'a':a,
                                               'Qbb':Qbb,
                                               'Sigma_bb':Sigma_bb,
                                               'bb':bb,
                                               'bbtilde':bbtilde}            

dumped = json.dumps(optics, cls=NumpyEncoder)
with open('/Users/jkravz311/Desktop/EAP_optics.json', 'w') as fp:
    json.dump(dumped, fp)       

# for loading saved json w/ numpy arrays...
# with open('/Users/jkravz311/Desktop/EAP_optics.json', 'r') as f:
#     loaddata = json.load(f)
#     data = json.loads(loaddata)


