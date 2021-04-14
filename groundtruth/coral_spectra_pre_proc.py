#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pre-processing raw benthic spectra
Run cells individually

Cell 1: setup for tab delimited data
Cell 2: proc for tab delimited data
Cell 3: proc for csv delimited data

@author: jkravitz 2021
"""
#%% CELL 1: SETUP FOR TAB DELIMITED DATA

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import BSpline, LSQUnivariateSpline

# setup for all formats
path = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/corals/raw/Roelfsema_phinn/cook_islands_compiled_raw.csv'
outpath = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/corals/final/Roelfsema_2017_cook.csv'
# define knot positions
# minus first and last (400,800/900)
l = np.arange(400,801,1)
kp = [410,420,425,430,440,460,480,500,515,530,545,560,575,590,600,610,
      620,630,640,650,665,680,695,710,725,740,760,780]

# extra setup if tab delimited
# header = 82
# cite = 'Roelfsema & Phinn, 2013'
# loc = 'Heron reef, Australia'
# metaid = np.array(['Name','Comment','Location','Citation'])
# group1 = 'Event'
# group2 = 'Comment'
# quick check of data
# data = pd.read_csv(path,sep='\t',header=header)

#%% CELL 2: PROC FOR TAB DELIMITED DATA

def main (path,outpath,header,cite,loc,metaid,l,group1,group2,kp):
    
    dataout = pd.DataFrame()
    data = pd.read_csv(path,sep='\t',header=header)
    grouped = data.groupby([group1,group2])
    count = 0
    for k, group in grouped:
        print (k)
        lam = group['Lambda [nm]'].values
        ref = group['Refl tot'].values
        meta = [k[0],k[1],loc,cite]
        if len(np.where(ref<0)[0]) > 300:
            continue
        # clean and grid
        newref = clean(l,lam,ref)
        # smooth
        ysmooth = smooth(l, newref, kp)
        yfinal = np.concatenate([meta,ysmooth])
        dataout[count] = yfinal
        count += 1
    
    out = dataout.T
    cols = np.concatenate([metaid,l])
    out.columns = cols
    out.to_csv(outpath)
    return out


def interp_nans_zeros (y):
    # nans, x= nan_helper(y)
    # y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    ynew = np.where(y < 0, np.nan, y)
    ynew = np.where(ynew == 0, np.nan, ynew)
    return np.isnan(ynew), lambda z: z.nonzero()[0]


# clean and interpolate to 400:900
def clean (l, lam, y):
    #ynew = np.interp(l, lam, y)
    nans, x = interp_nans_zeros(y)
    y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    ynew = np.interp(l,lam,y)
    # idx = np.nonzero(data)
    # x = np.arange(len(data))
    # interp = interp1d(x[idx],data[idx],fill_value="extrapolate")
    # ynew = interp(x)

    return ynew  

def smooth (l, ref, kp):
        
    # spline
    s1 = LSQUnivariateSpline(l, ref, kp) # spline fit
    kn = s1.get_knots() # knots
    kn = 3*[kn[0]] + list(kn) + 3*[kn[-1]]
    c = s1.get_coeffs()
    s2 = BSpline(kn,c,3)
    ysmooth = s2(l)
        
    return ysmooth    


if __name__ == '__main__':
    out = main(path,outpath,header,cite,loc,metaid,l,group1,group2,kp)


# quick plot for check
spectra = out.iloc[:,4:]
spectra = spectra.astype(float)
spectra.T.plot(legend=False)

#%% CEll 3: For non tab delimeted data: just needs smoothing and interping

dataout = pd.DataFrame()
data = pd.read_csv(path)
count = 0
for k in data.index:
    info = data.loc[k,:]
    meta = info.iloc[0:4]
    ref = info.iloc[4:].astype(float)
    lam = ref.index.values.astype(float)
    # clean and grid
    newref = clean(l,lam,ref)
    # smooth
    ysmooth = smooth(l, newref, kp)
    yfinal = np.concatenate([meta,ysmooth])
    dataout[count] = yfinal
    count += 1
        
out = dataout.T
cols = np.concatenate([meta.index,l])
out.columns = cols    
out.to_csv(outpath)

# quick plot for check
spectra = out.iloc[:,4:]
spectra = spectra.astype(float)
fig, ax = plt.subplots()
spectra.T.plot(ax=ax,legend=False)
ax.set_ylim(0,4)

        