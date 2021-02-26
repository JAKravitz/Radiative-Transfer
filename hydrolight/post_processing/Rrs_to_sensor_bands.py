##
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from Rrs_to_SRF_funs import *
from sample_info import sample_info
from band_indices import *
import pickle
# ------------- #
db = 'cyano'
# ------------- #

flPath = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN//{0}/fluor/csv/'.format(db)
noFlPath = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/noFluor/csv/'.format(db)

# flPath = '/Volumes/500gb/hy_data/NN_fluor/csv/'
# noFlPath = '/Volumes/500gb/hy_data/NN_nofluor/csv/'

hicowv = pd.read_csv('/Users/jkravz311/Desktop/PhD/srf_resample/hicowv.csv',header=None)
dirlist = os.listdir(flPath)
dirlist2 = []
for f in dirlist:
    dirlist2.append(os.path.splitext(f)[0])
#dirlist2= dirlist2[1:]
wv = np.arange(400,901,1)
wvhico = hicowv.values[0]
l8wv = np.array([443,482,561,665,865])
s2wv = np.array([443,490,560,665,705,740,783,842,865])
modiswv = np.array([412,443,469,488,531,551,555,645,667,678,748,859,869])
meriswv = np.array([412.5,442.5,490,510,560,620,665,681.25,708.75,753.75,760.625,778.75,865,885,900])
s3wv = np.array([400,412.5,442.5,490,510,560,620,665,673.75,681.25,708.75,753.75,761.25,764.375,767.75,778.75,865,885])
l8srfpath = '/Users/jkravz311/Desktop/PhD/srf_resample/l8_srfs.csv'
modissrfpath = '/Users/jkravz311/Desktop/PhD/srf_resample/modis_srfs.csv'
merissrfpath = '/Users/jkravz311/Desktop/PhD/srf_resample/meris_srfs2.csv'
olcisrfpath = '/Users/jkravz311/Desktop/PhD/srf_resample/olci_srfs.csv'
olcisrfoptpath = '/Users/jkravz311/Desktop/PhD/srf_resample/olci_opt_srfs.csv'
s2srfpath = '/Users/jkravz311/Desktop/PhD/srf_resample/s2a_srfs.csv'
cysiop = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/cyano_siops.csv')

##
print 'Resampling SRFs...'
hicorrs = pd.DataFrame(columns=wvhico,index=dirlist2)
l8rrs = pd.DataFrame(columns=dirlist2)
modisrrs = pd.DataFrame(columns=dirlist2)
merisrrs = pd.DataFrame(columns=dirlist2)
s3rrs = pd.DataFrame(columns=dirlist2)
s2rrs = pd.DataFrame(columns=dirlist2)
hyperrs = pd.DataFrame(columns=dirlist2)

for f in dirlist2:

    info = f.split('_')

    # interp to 1nm resolution
    data = pd.read_csv(os.path.join(flPath,f+'.csv'),index_col=0)
    rrs_intp = np.interp(wv,data.index,data.Rrs)
    data2 = pd.DataFrame(index=wv,data=rrs_intp)
    hyperrs[f] = rrs_intp

    # to hico spec response functions
    for k in wvhico:
        if k < 745:
            b1 = round(k - 5)
            b2 = round(k + 5)
            mn = np.mean(data2.loc[b1:b2])
            hicorrs.at[f,k] = mn.values[0]
        else:
            b1 = round(k - 10)
            b2 = round(k + 10)
            mn = np.mean(data2.loc[b1:b2])
            hicorrs.at[f,k] = mn.values[0]
    # l8 spec functions
    rrs = resample_l8(data2, l8srfpath)
    l8rrs[f] = rrs.values[:,0]
    # modis
    rrs = resample_modis(data2, modissrfpath)
    modisrrs[f] = rrs.values[:,0]
    # meris
    rrs = resample_meris(data2, merissrfpath)
    merisrrs[f] = rrs.values[:,0]
    # olci
    rrs = resample_OLCI(data2,olcisrfpath)
    s3rrs[f] = rrs.values[:,0]
    # s2
    rrs = resample_s2(data2,s2srfpath)
    s2rrs[f] = rrs.values[:,0]

#
l8rrs.index = l8wv
modisrrs.index = modiswv
merisrrs.index = meriswv
s3rrs.index = s3wv
s2rrs.index = s2wv
hyperrs.index= wv
#
l8rrs = l8rrs.T
modisrrs = modisrrs.T
merisrrs = merisrrs.T
s3rrs = s3rrs.T
s2rrs = s2rrs.T
hyperrs = hyperrs.T
#
sensorDict = {'l8':l8rrs,'modis':modisrrs,'meris':merisrrs,'s3':s3rrs,'s2':s2rrs,'hico':hicorrs,'hyper':hyperrs}

f = open('/Users/jkravz311/Desktop/dict1.pkl','wb')
pickle.dump(sensorDict,f)
f.close()

## SAMPLE INFO

fin = open('/Users/jkravz311/Desktop/dict1.pkl','rb')
sensorDict = pickle.load(fin)

print 'Combining sample info...'

sensorDict2 = sample_info(flPath,noFlPath,sensorDict)

f = open('/Users/jkravz311/Desktop/dict2.pkl','wb')
pickle.dump(sensorDict2,f)
f.close()

## SENSOR INDICES
# import pickle
# fin = open('/Users/jkravz311/Desktop/dict2.pkl','rb')
# sensorDict2 = pickle.load(fin)

print 'Calculating band indices...'
sensorDict3 = band_indices(sensorDict2)

f = open('/Users/jkravz311/Desktop/dict3.pkl','wb')
pickle.dump(sensorDict3,f)
f.close()

##
sensorDict3['meris'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/MERISrrs.csv'.format(db))
sensorDict3['modis'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/MODISrrs.csv'.format(db))
sensorDict3['l8'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/l8rrs.csv'.format(db))
sensorDict3['hico'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/HICOrrs.csv'.format(db))
sensorDict3['s3'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/s3rrs.csv'.format(db))
sensorDict3['s2'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/s2rrs.csv'.format(db))
sensorDict3['hyper'].to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/HYPErrs.csv'.format(db))


##
# merisrrs.to_csv('/Volumes/500gb/{0}/merisrrs.csv'.format(db))
# modisrrs.to_csv('/Volumes/500gb/{0}/modisrrs.csv'.format(db))
# l8rrs.to_csv('/Volumes/500gb/{0}/l8rrs.csv'.format(db))
# hicorrs.to_csv('/Volumes/500gb/{0}/hicorrs.csv'.format(db))
# s3rrs.to_csv('/Volumes/500gb/{0}/s3rrs.csv'.format(db))
# s3rrs_opt.to_csv('/Volumes/500gb/{0}/s3optrrs.csv'.format(db))
# s2rrs.to_csv('/Volumes/500gb/{0}/s2rrs.csv'.format(db))
# #
# s3rrs_idx_final.to_csv('/Volumes/500gb/{0}/s3rrs_idx.csv'.format(db))
# s3rrs_opt_idx_final.to_csv('/Volumes/500gb/{0}/s3optrrs_idx.csv'.format(db))
# s2rrs_idx_final.to_csv('/Volumes/500gb/{0}/s2rrs_idx.csv'.format(db))
# merisrrs_idx_final.to_csv('/Volumes/500gb/{0}/merisrrs_idx.csv'.format(db))
# modisrrs_idx_final.to_csv('/Volumes/500gb/{0}/modisrrs_idx.csv'.format(db))

print 'FINISHED!'