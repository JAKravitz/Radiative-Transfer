##
import os
import re
import pandas as pd
import numpy as np

# db = 'cyano'
# fl = 'fluor'

# path = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/{1}/txt/'.format(db,fl)
# outpath = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/{0}/{1}/csv/'.format(db,fl)

path = '/Volumes/500gb/hy_data/NN_fluor/txt/'
outpath = '/Volumes/500gb/hy_data/NN_fluor/csv/'

flist = os.listdir(path)
print 'Running ...'
for file in flist:

    if not file.startswith('M'):
        continue
    name = os.path.splitext(file)[0]
    fname = os.path.join(path,file)

    wl = []
    rrs = []
    rrsfile = pd.DataFrame()
    eu = []
    ed = []
    eo = []
    lu = []
    with open(fname, 'r+') as f:
        for i,line in enumerate(f):

            if i == 0:
                info = line
                title = re.findall(r'I_.*txt', info)

            if 384 < i < 470:
                wl.append(line.split(' ')[1])
                rrs.append(line.split(' ')[4])

    #
    newname = os.path.splitext(title[0])[0]
    newname = newname[2:]

    wl = np.array(wl,dtype=float)
    rrs = np.array(rrs,dtype=float)
    rrsfile['Rrs'] = rrs
    rrsfile.index = wl
    rrsfile.to_csv(outpath+'{0}.csv'.format(newname))
    #

print 'DONE!'
##


