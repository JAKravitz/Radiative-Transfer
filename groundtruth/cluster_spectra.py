'''
Clustering spectra using beta spline coefficients
Elbow method for defining optimal # k for Kmeans clustering

Kravitz, 2021
'''

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skfda.representation.basis as basis
from scipy.interpolate import BSpline, LSQUnivariateSpline
from skfda.exploratory.visualization import Boxplot
from skfda.representation.grid import FDataGrid
from matplotlib.colors import LinearSegmentedColormap
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

### INPUTS ###
# path to compiled aster data
#datapath = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/jpl_aster_data.csv'

##############


def standardize_spectra (data):
    st_data = pd.DataFrame()
    for k in data.index:
        spec = data.loc[k,:]
        stspec = spec / np.trapz(y=spec.values,x=spec.index.astype(int),axis=0)
        st_data[k] = stspec
    st_data = st_data.T
    return st_data
    
def beta_spline (data,t2):
    # data should be standardized spectra
    coefs = pd.DataFrame()
    for k in data.index:
        spec = data.loc[k,:]
        x = spec.index.values.astype(int) # x-data
        y = spec.values # y-data
        s1 = LSQUnivariateSpline(x, y, t2) # spline fit
        kn = s1.get_knots() # knots
        kn = 3*[kn[0]] + list(kn) + 3*[kn[-1]]
        c = s1.get_coeffs()
        s2 = BSpline(kn,c,3)
        c2 = list(s2.c)
        coefs[k] = c2
    coefs = coefs.T
    return coefs

def spline_coefs (data):
    
    # load data and separate spectra
    #data = pd.read_csv(datapath)
    spectra = data.filter(regex='\d')
    meta = data.filter(regex='\D')

    # define knots
    knots = [400,430,460,500,515,530,545,560,575,590,605,620,635,650,665,680,695,710,725,740,755,800,850,900]
    # define knot positions
    t2 = [430,460,500,515,530,545,560,575,590,605,620,635,650,665,680,695,710,725,740,755,800,850]

    # plot basis functions for verification
    # bss = basis.BSpline([400,900],knots=knots)
    # fig, ax = plt.subplots()
    # bss.plot(chart=ax)
    # ax.set_xlabel('Wavelength (nm)')
    # fig.savefig('/Users/jkravz311/Desktop/bspline_funs.png',bbox_inches='tight',dpi=300)

    # standardize original spectra
    st_data = standardize_spectra(spectra)

    # convert spectra to basis functions
    coefs = beta_spline(st_data,t2)
    
    return coefs, spectra, meta


def elbow (data):
    # elbow method to define # clusters
    # requires visual inspection of plot
    data = data.values(data)
    K = range(1,20)
    dist = []
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(data)
        dist.append(km.inertia_)
    fig, ax = plt.subplots()
    ax.plot(K,dist,'ko-')


def kmeans_clust (data,k):
    # k means clustering
    kmeans = KMeans(n_clusters=k)
    kmeans = kmeans.fit(data)
    labels = kmeans.predict(data)
    #centroids = kmeans.cluster_centers_
    return labels


#%%
# spectra['labels'] = labels
# clusters = spectra.groupby('labels')
# clustDict = {}
# for key, group in clusters:
#     cluster = group.iloc[:,:-1]
#     clustDict[key] = cluster
#     cluster.T.plot(legend=False)
#     print (key)


