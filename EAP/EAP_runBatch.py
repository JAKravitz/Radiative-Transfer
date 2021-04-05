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
import pickle
from EAP_proc import EAP


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

# paths to imaginary and real refractive index
start = 'P. calceolata' # 0 to start from beginning, else phyto
mf = '/Users/jkravz311/git_projects/Radiative-Transfer/EAP/501nm_extended_e1701000.mat'
astarpath = '/Users/jkravz311/git_projects/Radiative-Transfer/EAP/data/in_vivo_phyto_abs.csv'
batchinfo = pd.read_csv('/Users/jkravz311/git_projects/Radiative-Transfer/EAP/data/EAP_batch_V1.csv', index_col=0)

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
          'nshell1': np.arange(1.06,1.16,.01),
          'nshell2': np.arange(1.11, 1.22,.01),
          'nshell3': np.arange(1.1, 1.19,.01),
          'ncore': np.arange(1.014, 1.04, 0.005)}

#%%

## if json file...
# class NumpyEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, np.ndarray):
#             return obj.tolist()
#         return json.JSONEncoder.default(self, obj)


# dumped = json.dumps(optics, cls=NumpyEncoder)
# with open('/Users/jkravz311/Desktop/EAP_optics.json', 'w') as fp:
#     json.dump(dumped, fp) 

def pandafy (array, Deff):
    out = pd.DataFrame(array, index=Deff)
    return out

# define where to start in batch list
if start == 0:
    with open('/Users/jkravz311/Desktop/EAP_optics.p', 'wb') as fp:
        pickle.dump(optics, fp) 
else:
    batchinfo = batchinfo.loc[start:,:]

# loop through phyto batch list
for i,phyto in enumerate(batchinfo.index):
    
    # with open('/Users/jkravz311/Desktop/EAP_optics.json', 'r+') as fp:
    #     loaddata = json.load(fp)
    #     data = json.loads(loaddata)

    with open('/Users/jkravz311/Desktop/EAP_optics.p', 'rb') as fp:
        data = pickle.load(fp)
        
        print ('####### i: {} - phyto: {} #######'.format(i,phyto))
        
        # get sample info
        info = batchinfo.loc[phyto,:]
        clss = info.Class
        VsF = np.random.choice(params[info.Vs], 3, replace=False)
        CiF = np.random.choice(params[info.Ci], 3, replace=False)
        nshellF = np.random.choice(params[info.nshell], 3, replace=False)
        ncore = np.random.choice(params['ncore'], 1)
        Deff = np.arange(info.Dmin, info.Dmax, 1)
        
        # add new phtyo to dictionary
        data[clss][phyto] = {}
        
        # loop through parameters
        for Vs in VsF:
            for ci in CiF:
                for n in nshellF:
                    
                    # run name format: 'Vs_Ci_nshell'
                    rname = '{:.2f}_{:.2f}_{:.2f}'.format(Vs, ci/1e6, n)
                    data[clss][phyto][rname] = {}
                    
                    # EAP run
                    # standard
                    print ('####### i: {} - phyto: {} #######'.format(i,phyto))
                    print ('------ {} ------'.format(rname))
                    Qc, Sigma_c, c, Qb, Sigma_b, b, Qa, Sigma_a, a, Qbb, Sigma_bb, bb, bbtilde = EAP(phyto, mf, astarpath, Vs, ci, Deff, n, ncore)
                    
                    # empty dict for current run
                    rname_data = {'Qc': Qc,
                                  'Sigma_c': Sigma_c,
                                  'c': c,
                                  'Qb': Qb,
                                  'Sigma_b': Sigma_b,
                                  'b': b,
                                  'Qa': Qa,
                                  'Sigma_a': Sigma_a,
                                  'a': a,
                                  'Qbb': Qbb,
                                  'Sigma_bb': Sigma_bb,
                                  'bb': bb,
                                  'bbtilde': bbtilde}  
                    
                    # pandafy params so Deff is index
                    for param in ['Qc','Sigma_c','c','Qb','Sigma_b','b','Qa','Sigma_a','a','Qbb','Sigma_bb','bb','bbtilde']:
                        rname_data[param] = pandafy(rname_data[param], Deff)
        
                    # save current run to pickle dict
                    data[clss][phyto][rname] = rname_data 
           
    # dumped = json.dumps(data, cls=NumpyEncoder)
    # json.dump(dumped,fp)
    with open('/Users/jkravz311/Desktop/EAP_optics.p', 'wb') as fp:
        pickle.dump(data,fp)

# for loading saved json w/ numpy arrays...
# with open('/Users/jkravz311/Desktop/EAP_optics.json', 'r') as f:
#     loaddata = json.load(f)
#     data = json.loads(loaddata)

#%% test
import pickle

with open('/Users/jkravz311/Desktop/EAP_optics.p', 'rb') as fp:
    foo = pickle.load(fp)



