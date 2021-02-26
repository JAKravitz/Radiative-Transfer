##
from __future__ import division
import os
import pandas as pd
import numpy as np
from lognorm_params import *

''' SET HYDROLIGHT PARAMETERS '''

sname_title = 'Case1'

# chl
sigma,scale = lognorm_params(1,5)
chlaData = lognorm_random(sigma, scale, 15000)

# admix
admix = 0

# cnap
sigma, scale = lognorm_params(.05, .3)
cnapData = lognorm_random(sigma,scale, 10000)
cnapData = cnapData[cnapData<=2]

# dino diameter
dinoD = np.arange(5,25,1)

# depth
depth = np.array([0,1,2,3,4,5,6,7,8,9])
depth = depth.reshape(-1, 1)
# solar zenith
thetaData = np.arange(5,70,5)

# paths
sioppath = '/Users/jkravz311/.wine/drive_c/HE52_data/SIOPs/'
iopzpath = '/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/iopz_files/'.format(sname_title)
nnpath = '/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/NN_files/'.format(sname_title)
if not os.path.exists(iopzpath):
    os.makedirs(iopzpath)
if not os.path.exists(nnpath):
    os.makedirs(nnpath)

# spectral parameters
R1 = np.arange(0,1.05,.05)
R2 = np.arange(0,1.05,.05)
R4 = np.arange(0,1.05,0.05)
P4 = np.arange(0.06,0.6,.06)

s_cdom = np.linspace(0.012,0.021,5)
astarnap = np.linspace(0.02,.3,7)
s_nap = np.linspace(0.007,.015,5)
bstarnap = np.linspace(0.5,1,5)
y2 = np.linspace(0.5,2,5)

wl = np.arange(400,901,1)
idx440 = np.where(wl==440)

##
sname = []
sname2 = []
for k in range(50):

    chl = round(np.random.choice(chlaData), 2)
    cnap = round(np.random.choice(cnapData), 3)
    dino = np.random.choice(dinoD)
    theta = np.random.choice(thetaData)

    sname.append(
                 str(chl) + '_' + str(admix) + '_' + str(dino) + '_' + str(cnap) + '_' + str(theta)
    )

for nii in sname:

    info = nii.split('_')
    chl = float(info[0])
    admix = float(info[1])
    dinoD = int(info[2])
    cnap = float(info[3])
    theta = int(info[4])

    # Get IOPs
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

    # pre - assign
    a_tot = []
    a_totphy = []
    c_tot = []
    c_totphy = []
    bb_tot = []
    bb_totphy = []
    chlProfile = []

    # chl profile
    chlProfile.append([chl] * 10)
    chlProfile = np.asarray(chlProfile).reshape(-1, 1)

    # spectral components
    # abs
    astar_nap = np.random.choice(astarnap,1)
    s_dm = np.random.choice(s_nap, 1)
    s_y = np.random.choice(s_cdom,1)
    # scatter
    bstar_nap = np.random.choice(bstarnap,1)
    y_2 = np.random.choice(y2,1)

    for i in range(10):

        # ABSORPTIONS
        # aphy
        aphy = (chlProfile[i] * (1-admix) * astard)
        aphy440 = aphy[idx440]

        # det/minerals absorption - NAP
        anap440 = cnap * astar_nap
        anap = anap440*np.exp(-s_dm*(wl-440))

        # glebstoff absorption
        p2 = 0.3 + (
                    (5.7 * np.random.choice(R2,1) * aphy440) /
                    (0.02 + aphy440)
                    )
        ag440 = p2 * aphy440
        ag = ag440 * np.exp(-s_y*(wl-440))

        # total absorption
        a_tot.append(aphy + anap + ag)
        a_totphy.append(aphy)

        # SCATTERING
        # bphy
        bphy = (chlProfile[i] * (1 - admix) * bstard)

        # bnap
        bnap550 = bstar_nap * cnap
        bnap = bnap550*(550/wl)**y_2

        # total attenuation (scatter + absorption)
        c_tot.append(aphy + anap + ag + bphy + bnap)
        c_totphy.append(aphy + bphy)

        # BACKSCATTERING
        # bbphy
        bbphy = (chlProfile[i] * (1 - admix) * bbstard)

        # bbnap
        bbnap = bnap * 0.0183

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

    ag440 = round(ag440,3)
    # WRITE TO .TXT
    nii = str(chl) + '_' + str(admix) + '_' + str(dinoD) + '_' + str(cnap) + '_' + \
            str(ag440) + '_' + str(theta)
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

    # IOPS ------------------------------------------------------

    aphy440 = aphy[40]
    bphy440 = bphy[40]
    bnap440 = bnap[40]
    bbphy440 = bbphy[40]
    bbnap440 = bbnap[40]

    nn_output = pd.Series()
    nn_output['Chl'] = chl
    nn_output['Admix'] = admix
    nn_output['PC_li'] = np.nan
    nn_output['PC_simis'] = np.nan
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
    sname_df.at[k, 'cy'] = np.nan

sname_df.to_csv('/Users/jkravz311/.wine/drive_c/HE52_data/NN/{0}/sname_{0}.csv'.format(sname_title))

print 'All Chlz, PCz, and IOPz files written.'

##
