##
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import cos
from astropy.modeling.models import Gaussian1D

#hyfilespath = '/Users/jkravz311/.wine/drive_c/HE52/output/Ecolight/excel/'
hyfilespath = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/cyano/noFluor/txt/'
# path to where new HE txt files will go with fluorescence added
newpath = '/Users/jkravz311/Desktop/PhD/cyano_sims/NN/cyano/fluor/txt/'

flist = os.listdir(hyfilespath)
wv = np.arange(400,900.5,.5)
wv2 = np.arange(400,901,1)

#----------
fqyList = [0.005,0.01,0.015,0.02]
QaList = [0.3,0.4,0.5,0.6]
#----------

##
def fl_calc(hyfilespath,fname,newpath,wv,wv2):

##
    #fname = 'M0.17_0.4_8_0.485_1.162_45_aer.txt'
    import pandas as pd
    import numpy as np
    from scipy import integrate

    path = hyfilespath + fname
    astarpath = '/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/'

    # # cyano
    # info = fname.split('_')
    # cy = info[7].split('.')[0]
    # cy = 'aer'

    # get wavelengths
    wl = []
    with open(path, 'r+') as f:
        for i, line in enumerate(f):
            if 22 < i < 108:
                info = line.split('  ')
                wl.append(info[0])
    wl2 = np.array(wl).astype(float)

    # get optical data required for fluorescence calculation from hydrolight output
    btot = pd.DataFrame(index=wv, columns=[0,1,2,3,4])
    atot = pd.DataFrame(index=wv, columns=[0,1,2,3,4])
    rrs = []

    with open(path, 'r+') as f:
        for i, line in enumerate(f):

            # run title
            if i == 0:
                info = line
                title = re.findall(r'I_.*txt', info)

            # Kpar
            if i == 7:
                info = line.split('  ')
                Kpar = float(info[4])

            # PAR
            if i == 15:
                info = line.split('  ')
                PAR = float(info[4])

            # Rrs
            if 374 < i < 460:
                rrs.append(line.split(' ')[4])

            # btot
            if 1782 < i < 1868:
                info = line.split('  ')
                for k, d in enumerate(btot.columns.values):
                    btot.at[float(info[0]), d] = info[k + 1]
            btot = btot.astype(float)

            # atot
            if 2046 < i < 2132:
                info = line.split('  ')
                for k, d in enumerate(atot.columns.values):
                    atot.at[float(info[0]), d] = info[k + 1]
            atot = atot.astype(float)

    # interpolate
    btot.interpolate(limit_direction='both',limit=8,inplace=True); btot = btot.loc[wv2,:]
    atot.interpolate(limit_direction='both',limit=8,inplace=True); atot = atot.loc[wv2,:]
    ctot = atot + btot

    rrs = np.array(rrs).astype(float)
    rrs = np.interp(wv2,wl2,rrs)

    # get astarphy dinos and cyanos
    info = os.path.splitext(title[0])[0]
    info = info[2:]
    info = info.split('_')
    chl = float(info[0])
    admix = float(info[1])
    dsize = info[2]
    cy = info[6]

    astar = pd.read_csv(astarpath + 'astarphy_dinos/astarphy_dino_{0}.csv'.format(dsize), index_col=0, header=None)
    astar_cy = pd.read_csv(astarpath + 'cyano_siops.csv',index_col=0)
    astar_cy = astar_cy.loc[:,'a_{0}'.format(cy)]

    # depth integration over optical depth

    # astarphy
    astarcy = np.interp(wv2,astar_cy.index, astar_cy.values)
    astarcy = pd.Series(astarcy * (chl * admix), index=wv2)
    astarcy = pd.concat([astarcy] * atot.shape[1], axis=1)
    astarcy.columns = atot.columns

    astarphy = np.interp(wv2, astar.index, astar.iloc[:, 0].values)
    astarphy = pd.Series(astarphy * (chl*(1-admix)), index=wv2) #

    astarphy = pd.concat([astarphy] * atot.shape[1], axis=1)
    astarphy.columns = atot.columns

    # zmax = round(zmax685,0)

    astarphy2 = astarphy.iloc[0:300,0]
    astarcy2 = astarcy.iloc[0:300,0]

    ctot685_g1 = ctot.loc[685,:]
    ctot685_g2 = ctot.loc[730,:]

    ctot6852_g1 = ctot685_g1.loc[0]
    ctot6852_g2 = ctot685_g2.loc[0]

    # attenuated
    Ktot2_g1 = Kpar + ctot6852_g1
    Ktot2_g2 = Kpar + ctot6852_g2

    # absorbed
    absdphy = integrate.simps(astarphy2.values, astarphy2.index) * PAR + \
              ((integrate.simps(astarcy2.values, astarcy2.index) * .15) * PAR)
    absd_g1 = absdphy / Ktot2_g1
    absd_g2 = absdphy / Ktot2_g2

    # calculation of fl amplitude
    fqy = np.random.choice(fqyList)
    Qa = np.random.choice(QaList)
    lf_g1 = (.54 * (1 / (4 * np.pi)) * (fqy/25) * Qa) * absd_g1
    lf_g2 = (.54 * (1 / (4 * np.pi)) * (fqy/50) * Qa) * absd_g2
    flamp_g1 = (lf_g1 / 1.23) * .001
    flamp_g2 = (lf_g2 / 1.23) * .001

    # gaussian function - bimodal
    g1 = Gaussian1D(amplitude=flamp_g1, mean=685, stddev=25/2.355) + Gaussian1D(amplitude=flamp_g2, mean=730, stddev=50/2.355)
    fl = g1(wv2)
    fl2 = np.interp(wl2,wv2,fl)

    # Preview Plot
    # fig, ax = plt.subplots()
    # ax.plot(wv2,fl, color='green', label='fl')
    # ax.plot(wv2,rrs,color='blue',label='rrs')
    # ax.plot(wv2,rrs + fl, color='red',label='total')
    # ax.set_xlim(400,800)
    # ax.legend()
    # ax.set_ylabel('Rrs')
    # ax.set_xlabel('Wavelength (nm)')
    # fig.savefig('/Users/jkravz311/Desktop/fl.png',bbox_inches='tight',dpi=125)

##
    # BIMODAL DISTRIBUTION OF CHL FL

    # def gauss(x, mu, sigma, A):
    #     return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)
    #
    # def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    #     return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)
    #
    # fl1 = flamp_g1
    # fl2 = flamp_g2
    # s1 = 25 / 2.355
    # s2 = 50 / 2.355
    # params = [685, s1, fl1, 730, s2, fl2]
    #
    # fig, ax = plt.subplots()
    # ax.plot(wv2, bimodal(wv2, *params), color='green', label='fl')
    # ax.plot(wv2,rrs, color='blue', label = 'rrs')
    # ax.plot(wv2,rrs + bimodal(wv2,*params), color = 'red', label='total')
    # ax.set_xlim(400,800)
    # ax.legend()
    # ax.set_ylabel('Rrs')
    # ax.set_xlabel('Wavelength (nm)')
    # fig.savefig('/Users/jkravz311/Desktop/fl.png',bbox_inches='tight',dpi=125)

    #
    newfname = 'M' + '_'.join(info) + '_{0}.txt'.format(fqy)
    #

    def replace_all(text, dic):
        for i, j in dic.iteritems():
            text = text.replace(i, j)
        return text

    with open(newpath+newfname, 'w') as newf:
        with open(path, 'r') as f:
            count = 0
            for i, line in enumerate(f):
                if 374 < i < 460:
                    #print line.split('  ')
                    rrs = float(line.split('  ')[1])
                    ed = float(line.split('  ')[2])
                    lw = float(line.split('  ')[3])
                    lu = float(line.split('  ')[4])
                    le = lu-lw

                    newrrs = ' ' + str(rrs+fl2[count])
                    newlw = ' ' + str((rrs+fl2[count]) * ed)
                    newlu = ' ' + str(((rrs + fl2[count]) * ed) +le) + '\r\n'

                    d = {line.split('  ')[1]: newrrs,
                         line.split('  ')[3]: newlw,
                         line.split('  ')[4]: newlu}

                    # newf.write(line.replace(line.split('  ')[1], newrrs))
                    # newf.write(line.replace(line.split('  ')[3], newlw))
                    # newf.write(line.replace(line.split('  ')[4], newlu))

                    newline = replace_all(line,d)

                    newf.write(line.replace(line,newline))

                    count = count+1
                else:
                    newf.write(line)

##

print 'running...'

for fname in os.listdir(hyfilespath):

    fl_calc(hyfilespath,fname,newpath,wv,wv2)

print 'DONE!'



