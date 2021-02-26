##
from __future__ import division
import os
import pandas as pd
import numpy as np
from lognorm_params import *

''' SET HYDROLIGHT PARAMETERS '''

sname_title = 'CDOM'

# chl
sigma,scale = lognorm_params(10,50)
chlaData = lognorm_random(sigma, scale, 40000)
chlaData = chlaData[chlaData<=150]

# admix
admix0 = 0
admix1 = np.arange(0,.8,.1)
admix2 = np.array([0,.1,.2,.3,])

# dino diameter
dino1 = np.arange(5,20,1)
dino2 = np.arange(15,35,1)
dino3 = np.arange(30,45,1)

# cnap
sigma, scale = lognorm_params(1, 4)
cnapData = lognorm_random(sigma,scale, 10000)
cnapData = cnapData[cnapData<=20]

# cdom
sigma, scale = lognorm_params(5, 10)
acdomData = lognorm_random(sigma,scale, 10000)
acdomData = acdomData[acdomData<=25]

# depth
depth = np.array([0,1,2,3,4,5])
depth = depth.reshape(-1, 1)

# solar zenith
thetaData = np.arange(5,70,5)

sname = []
sname2 = []
for k in range(20000):

    # run info
    chl = round(np.random.choice(chlaData), 2)
    cnap = round(np.random.choice(cnapData), 3)
    cdom = round(np.random.choice(acdomData), 3)
    theta = np.random.choice(thetaData)


    if 1 <= cdom < 5:

        admix = np.random.choice(admix1)
        if 0 <= chl < 20:
            dino = np.random.choice(dino1)
        elif 20 <= chl < 50:
            dino = np.random.choice(dino2)
        elif chl >= 50:
            dino = np.random.choice(dino3)

        sname.append(
            str(chl) + '_' + str(admix) + '_' + str(dino) + '_' + str(cnap) + '_' + str(cdom) + '_' + str(theta)
        )

    elif 5 <= cdom < 10:

        admix = np.random.choice(admix2)
        if 0 <= chl < 20:
            dino = np.random.choice(dino1)
        elif 20 <= chl < 50:
            dino = np.random.choice(dino2)
        elif chl >= 50:
            dino = np.random.choice(dino3)

        sname.append(
            str(chl) + '_' + str(admix) + '_' + str(dino) + '_' + str(cnap) + '_' + str(cdom) + '_' + str(theta)
        )

    elif cdom >= 10:

        admix = np.random.choice(admix0)
        if 0 <= chl < 20:
            dino = np.random.choice(dino1)
        elif 20 <= chl < 50:
            dino = np.random.choice(dino2)
        elif chl >= 50:
            dino = np.random.choice(dino3)

        sname.append(
            str(chl) + '_' + str(admix) + '_' + str(dino) + '_' + str(cnap) + '_' + str(cdom) + '_' + str(theta)
        )

##
''' WRITE IOP FILES'''

print 'Writing IOP, Chlz, and PCz files...'

sioppath = '/Users/jkravz311/.wine/drive_c/HE52_data/SIOPs/'
iopzpath = '/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/iopz_files/'.format(sname_title)
nnpath = '/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/NN_files/'.format(sname_title)
if not os.path.exists(iopzpath):
    os.makedirs(iopzpath)
if not os.path.exists(nnpath):
    os.makedirs(nnpath)

# spectral parameters
s_cdom = np.linspace(0.012,0.021,5)
astarnap = np.linspace(0.02,.3,7)
s_nap = np.linspace(0.007,.015,5)
bstarnap = np.linspace(0.5,1,5)
y2 = np.linspace(0.5,2,5)

wl = np.arange(400,901,1)
idx440 = np.where(wl==440)

cyList = ['aer','ana','aph','nod']
for nii in sname:

    info = nii.split('_')
    chl = float(info[0])
    admix = float(info[1])
    dinoD = int(info[2])
    cnap = float(info[3])
    ag440 = float(info[4])
    theta = int(info[5])

    # get SIOPs
    # dinos a
    astard = pd.read_csv(sioppath+'dino_a.csv',index_col=0)
    astard = astard.loc[dinoD,:]
    astard = np.interp(wl,astard.index.values.astype(int),astard.values)
    # dinos b
    bstard = pd.read_csv(sioppath+'dino_b.csv',index_col=0)
    bstard = bstard.loc[dinoD,:]
    bstard = np.interp(wl,bstard.index.values.astype(int),bstard.values)
    # dinos bb
    bbstard = pd.read_csv(sioppath+'dino_bb.csv',index_col=0)
    bbstard = bbstard.loc[dinoD,:]
    bbstard = np.interp(wl,bbstard.index.values.astype(int),bbstard.values)
    # cyanos
    cysiop = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/cyano_siops.csv',index_col=0)
    cy = np.random.choice(cyList, 1)[0]
    astarcy = np.interp(wl, cysiop.index, cysiop.loc[:,'a_{0}'.format(cy)])
    bstarcy = np.interp(wl, cysiop.index, cysiop.loc[:,'b_{0}'.format(cy)])
    bbstarcy = np.interp(wl, cysiop.index, cysiop.loc[:,'bb_{0}'.format(cy)])

    # pre - assign
    a_tot = []
    a_totphy = []
    c_tot = []
    c_totphy = []
    bb_tot = []
    bb_totphy = []
    chlProfile = []

    # chl profile
    chlProfile.append([chl] * 6)
    chlProfile = np.asarray(chlProfile).reshape(-1, 1)

    # spectral components
    # abs
    astar_nap = np.random.choice(astarnap,1)
    s_dm = np.random.choice(s_nap, 1)
    s_y = np.random.choice(s_cdom,1)
    # scatter
    bstar_nap = np.random.choice(bstarnap,1)
    y_2 = np.random.choice(y2,1)

    for i in range(6):

        # ABSORPTIONS
        # aphy
        aphy = (chlProfile[i] * admix * astarcy) + (chlProfile[i] * (1-admix) * astard)
        # det/minerals absorption - NAP
        anap440 = cnap * astar_nap
        anap = anap440*np.exp(-s_dm*(wl-440))
        # glebstoff absorption
        ag = ag440*np.exp(-s_y*(wl-440))
        # total absorption
        a_tot.append(aphy + anap + ag)
        a_totphy.append(aphy)

        # SCATTERING
        # bphy
        bphy = (chlProfile[i] * admix * bstarcy) + (chlProfile[i] * (1 - admix) * bstard)
        # bnap
        bnap550 = bstar_nap * cnap
        bnap = bnap550*(550/wl)**y_2
        # total attenuation (scatter + absorption)
        c_tot.append(aphy + anap + ag + bphy + bnap)
        c_totphy.append(aphy + bphy)

        # BACKSCATTERING
        # bbphy
        bbphy = (chlProfile[i] * admix * bbstarcy) + (chlProfile[i] * (1 - admix) * bbstard)
        # bbnap
        bbnap = bnap * 0.02 # slightly higher than petzold (0.0183) for stronger mineral component (Gilerson07)
        # total backscatter
        bb_tot.append(bbphy + bbnap)
        bb_totphy.append(bbphy)

    # PREPARE FOR WRITING TO .TXT FILE IN HYDROLIGHT FORMAT
    l = len(wl)
    line11 = np.hstack((l,wl)).reshape(1,-1)
    # AC data
    foot = np.hstack((-1,[0]*(l*2)))
    # only phy
    acProfile_phy = np.hstack((depth,a_totphy,c_totphy))
    acProfile_phy = np.vstack((acProfile_phy,foot))
    # all components
    acProfile = np.hstack((depth,a_tot,c_tot))
    acProfile = np.vstack((acProfile,foot))

    # bb data
    foot = np.hstack((-1,[0]*(l)))
    # only phy
    bbProfile_phy = np.hstack((depth,bb_totphy))
    bbProfile_phy = np.vstack((bbProfile_phy,foot))
    # all components
    bbProfile = np.hstack((depth,bb_tot))
    bbProfile = np.vstack((bbProfile,foot))


    # WRITE TO .TXT
    nii = str(chl) + '_' + str(admix) + '_' + str(dinoD) + '_' + str(cnap) + '_' + \
            str(ag440) + '_' + str(theta) + '_{0}'.format(cy)
    sname2.append(nii)
    # acData
    fname = 'acData' + nii + '.txt'
    #fname_phy = 'acData' + nii + '_' + 'phy' + '.txt'
    header = '\n'.join(['AC DATA PROFILE',
                        'Format:chlz_admix_dinoD_cnap_acdom_fqy_apc_components',
                        'Does not include water!',
                        '#','#','#','#','#','#','#\n'])
    with open(os.path.join(iopzpath,fname),'w') as f:
        f.writelines(header)
        np.savetxt(f,line11,fmt='%d',delimiter='\t')
        np.savetxt(f,acProfile,delimiter='\t')

    # bbData
    fname = 'bbData' + nii + '.txt'
    #fname_phy = 'bbData' + nii + '_' + 'phy' + '.txt'
    header = '\n'.join(['BB DATA PROFILE',
                        'Format:chlz_admix_dinoD_cnap_acdom_fqy_apc_components',
                        'Does not include water!',
                        '#', '#', '#', '#', '#', '#', '#\n'])
    with open(os.path.join(iopzpath,fname),'w') as f:
        f.writelines(header)
        np.savetxt(f, line11, fmt='%d', delimiter='\t')
        np.savetxt(f, bbProfile, delimiter='\t')

##
# ------------------------------
    ''' WRITE NN OUTPUTS TO FILE '''
    # ------------------------------

    # PC ------------------------------------------------------

    from scipy import interpolate

    c_coef = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/apc_C_coeficients_li_2015.csv',
                         index_col=0)
    c1f = interpolate.interp1d(c_coef.index.values, c_coef.c1.values, fill_value='extrapolate')
    c1new = c1f(wl)

    c2f = interpolate.interp1d(c_coef.index.values, c_coef.c2.values, fill_value='extrapolate')
    c2new = c2f(wl)

    pc_Li = []
    pc_simis = []
    if admix == 0:
        pc_li = 0
        pc_sim = 0
    else:
        anw665 = a_tot[0][265]
        anw620 = a_tot[0][220]
        aphy_pc = 1.1872 * c1new * anw665 + c2new
        apc620 = anw620 - ag[220] - aphy_pc[220]
        astarpc = 0.008 * admix ** -.451  # from Matthews 2013
        pc_li = apc620 / astarpc

        apc620_sim = aphy[220] - 0.24 * aphy[265]
        pc_sim = apc620_sim / astarpc

    # for k in range(6):
    #     pc_Li.append(round(pc_li, 2))
    #     pc_simis.append(round(pc_sim, 2))

    # IOPS ------------------------------------------------------

    aphy440 = aphy[40]
    bphy440 = bphy[40]
    bnap440 = bnap[40]
    bbphy440 = bbphy[40]
    bbnap440 = bbnap[40]

    nn_output = pd.Series()
    nn_output['Chl'] = chl
    nn_output['Admix'] = admix
    nn_output['PC_li'] = round(pc_li,2)
    nn_output['PC_simis'] = round(pc_sim,2)
    nn_output['Cnap'] = cnap
    nn_output['aphy440'] = round(aphy440,4)
    nn_output['ag440'] = round(ag440,4)
    nn_output['anap440'] = round(anap440,4)
    nn_output['bphy440'] = round(bphy440,4)
    nn_output['bnap440'] = round(bnap440,4)
    nn_output['bbphy440'] = round(bbphy440,4)
    nn_output['bbnap440'] = round(bbnap440,4)

    NNname = 'NN_' + nii+'.csv'
    nn_output.to_csv(os.path.join(nnpath,NNname))


##

sname_df = pd.DataFrame(index=sname2, columns=['chla', 'admix', 'dinoD', 'cnap', 'cdom','theta','cy'])
for k in sname2:
    info = k.split('_')
    sname_df.at[k, 'chla'] = float(info[0])
    sname_df.at[k, 'admix'] = float(info[1])
    sname_df.at[k, 'dinoD'] = int(info[2])
    sname_df.at[k, 'cnap'] = float(info[3])
    sname_df.at[k, 'cdom'] = float(info[4])
    sname_df.at[k, 'theta'] = int(info[5])
    sname_df.at[k, 'cy'] = info[6]

sname_df.to_csv('/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/sname_{0}.csv'.format(sname_title))

print 'All Chlz, PCz, and IOPz files written.'

##
