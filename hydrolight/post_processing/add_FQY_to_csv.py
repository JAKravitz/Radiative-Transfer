##
import pandas as pd
import os

# flcsv = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/cyano/fluor/csv/'
# fltxt = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/cyano/fluor/txt/'
# data = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/cyano/s3rrs.csv',index_col=0)

flcsv = '/Volumes/500gb/hy_data/NN_fluor/csv/'
fltxt = '/Volumes/500gb/hy_data/NN_fluor/txt/'
data = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/others/s3rrs.csv',index_col=0)

##
print 'Running... '
txtlist = os.listdir(fltxt)

for k in txtlist:
    if not k.endswith('.txt'):
        continue
    info = k.split('_')
    name = info[:-1]
    fname = '_'.join(name)
    fqy = info[-1:][0]
    fqy = os.path.splitext(fqy)[0]
    #
    data.at[fname,'FQY'] = fqy


print 'DONE !!'

##
data1 = data.iloc[0:49858,:]
data2 = data.iloc[49858:,:]
data2 = data2.FQY
cols = data1.index.values
data2.index = cols
data1['FQY'] = data2

data1.to_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/NN/others/s3rrs.csv')