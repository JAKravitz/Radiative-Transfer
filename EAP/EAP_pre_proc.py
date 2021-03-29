#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pre-processing a*phy spectra for use as imaginary ref. index in EAP model

Remove nans and zeros from individual datasets
Resolve to standard lambda range
Smooth data using beta splines

@author: jkravitz 2021
"""
#%%
# smooth data using beta splines
# Remove nans and zeros from individual datasets
# Resolve to standard lambda range

from scipy.interpolate import griddata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline, LSQUnivariateSpline


path = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/phyto_optics/Clementson/C2.csv'
outpath = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/phyto_optics/Clementson/C2_smooth.csv'

def main (path,outpath):
    
    l = np.arange(.400, .901, .001)
    data = pd.read_csv(path, index_col=0)
    smoothspec, meta = smooth(data)
    smoothspec.index = data.index
    smoothspec.columns = l
    final = pd.concat([meta,smoothspec],axis=1)
    final.to_csv(outpath)

# def beta_spline_single (signal,kp):
    
#     x = signal
#     y = np.arange(400, 901, 1)
#     s1 = LSQUnivariateSpline(x, y, kp) # spline fit
#     kn = s1.get_knots() # knots
#     kn = 3*[kn[0]] + list(kn) + 3*[kn[-1]]
#     c = s1.get_coeffs()
#     s2 = BSpline(kn,c,3)
#     ysmooth = s2(x)
#     coefs = list(s2.c)
    
#     return coefs, ysmooth

def clean (data):
        
    # remove negs, nans, zeros
    s = np.where(data < 0, np.nan, data) # replace negatives w/ nan
    s = np.where(s == 0, np.nan, s) # replace zeros w/ nan
    s[np.isnan(s)] = min(s) # replaces all nans with minimum
    
    return s

def beta_spline_df (data,kp):
    
    coefs = pd.DataFrame()
    ysmooths = pd.DataFrame()
    for k in data.index:
        spec = data.loc[k,:]
        x = spec.index.values.astype(int) # x-data
        y = spec.values # y-data
        s1 = LSQUnivariateSpline(x, y, kp) # spline fit
        kn = s1.get_knots() # knots
        kn = 3*[kn[0]] + list(kn) + 3*[kn[-1]]
        c = s1.get_coeffs()
        s2 = BSpline(kn,c,3)
        ysmooth = s2(x)
        ysmooth = clean(ysmooth)
        c2 = list(s2.c)
        coefs[k] = c2
        ysmooths[k] = ysmooth
    coefs = coefs.T
    ysmooths = ysmooths.T
    
    return coefs, ysmooths


def smooth (data):
    
    newdata = []
    l = np.arange(400, 901, 1)
    
    # separate spectra and metadata
    spectra = data.filter(regex='\d')
    meta = data.filter(regex='[a-zA-Z]')
    wv = spectra.columns.values.astype(float)
    
    # define knot positions
    # minus first and last (400,900)
    kp = [410,420,425,430,440,460,480,500,515,530,545,560,575,590,600,610,
          620,630,640,650,665,680,695,710,725,740,760,780,800,850]
    
    for k in spectra.index:
        
        # resolve to standard 400-900nm @ 1nm res
        s = spectra.loc[k,:].values
        s = griddata(wv, s, l, 'linear') # interpolate 
        
        # clean
        s2 = clean(s)
        newdata.append(s2)
    
    # newdata to df
    data2 = np.array(newdata)
    data2 = pd.DataFrame(data2)
    data2.columns = l   

    # convert spectra to basis functions
    coefs, ysmooth = beta_spline_df(data2, kp)
        
    return ysmooth, meta


if __name__ == '__main__':
    main(path,outpath)




