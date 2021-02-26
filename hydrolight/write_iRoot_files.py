##
import os
import numpy as np
import pandas as pd
import shutil
import io

# PATHS:

# runlist = r'C:\HE52\run'
# IOPpath = r'C:\Users\User\Documents\Hy\data\iop_files'
# batchPath = r'C:\HE52\run\batch'
runlist = '/Users/jkravz311/.wine/drive_c/HE52/run/'
#IOPpath = '/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/iopz_files/'
batchPath = '/Users/jkravz311/.wine/drive_c/HE52/run/batch/'

# get filenames
fnames = pd.read_csv('/Users/jkravz311/.wine/drive_c/HE52_data/NN/sname.csv',index_col=0)
#fnames = pd.read_csv('/Users/jkravz311/Desktop/sname_2.csv',index_col=0)
sname = fnames.index.values

# Remove previous runlist for HE
if os.path.exists(os.path.join(runlist,'runlist.txt')):
    os.remove(os.path.join(runlist,'runlist.txt'))

# remove and replace batch
if os.path.exists(batchPath):
    shutil.rmtree(batchPath)
os.makedirs(batchPath)

##
''' --------- BEGIN WRITE NEW IROOT FILE --------- '''

#for fqy in [.02,.01,.005]:

for i,name in enumerate(sname):

    info = name.split('_')
    chl = float(info[0])
    admix = float(info[1])
    dinoD = float(info[2])
    cnap = float(info[3])
    cdom = float(info[4])


    ''' --------- SPECIFY LINES FOR IROOT FILE --------- '''

    # LINE 1 - RECORD 1: DEFAULT PARAMETERS
    # Specify quantum yield based on chl concentration
    if chl < 10:
        fqy = .01
        #line1 = np.array([0,400,700,0.01,488,0.00026,1,5.3]).reshape(1,-1)
    if chl >= 10:
        fqy = .005
        #line1 = np.array([0,400,700,0.005,488,0.00026,1,5.3]).reshape(1,-1)
    if chl >= 50:
        fqy = .003
        #line1 = np.array([0,400,700,0.003,488,0.00026,1,5.3]).reshape(1,-1)
    if chl >= 100:
        fqy = .002
        #line1 = np.array([0,400,700,0.002,488,0.00026,1,5.3]).reshape(1,-1)
    if chl >= 150:
        fqy = .001
        #line1 = np.array([0,400,700,0.001,488,0.00026,1,5.3]).reshape(1,-1)
    if chl >= 200:
        fqy = .0005
        #line1 = np.array([0, 400, 700, 0.0005, 488, 0.00026, 1, 5.3]).reshape(1, -1)
    if chl >= 300:
        fqy = .0001
        #line1 = np.array([0,400,700,0.0001,488,0.00026,1,5.3]).reshape(1,-1)

    line1 = np.array([0, 400, 700, fqy, 488, 0.00026, 1, 5.3]).reshape(1, -1)

    # create .txt filename
    fname0 = name #+ '_' + str(fqy)
    fname1 = 'I_' + fname0 + '_{0}.txt'.format(fqy)
    fname = os.path.join(batchPath,fname1)

    # LINE 2 - RECORD 2: RUN TITLE
    line2 = fname1+'\n'

    # LINE 3 - RECORD 3: ROOTNAME
    line3 = fname0 + '\n'

    # LINE 4 - RECORD 4A: OUTPUT OPTIONS
    line4 = np.array([-1, 0, 0, 1, 1, 1]).reshape(1, -1)

    # LINE 5 - RECORD 4B: MODEL OPTIONS
    line5 = np.array([3, 1, -1, 4, 0]).reshape(1, -1)

    # LINE 6 - RECORD 5A: NUMBER OF COMPONENTS
    line6 = np.array([2, 4]).reshape(1, -1)

    # LINE 7 - RECORD 5B: COMPONENT CONCENTRATIONS
    line7 = np.array([0, 0, 0, 0]).reshape(1, -1)

    # LINE 8-11 - RECORD 5C: SPECIFIC ABSORPTION PARAMETERS
    line8_11 = np.array([[0, 3, 440, 0.1, 0.014],  # water
                         [2, -666, 440, 0.1, 0.014],  # total
                         [0, -666, 440, 0.1, 0.017],  # cdom contribution
                         [2, 0, 440, 0.1, 0.014]]).reshape(4, -1)  # chl contribution

    # LINE 12-15 - RECORD 5D: SPECIFIC ABSORPTION DATA FILE NAMES
    line12_15 = '\n'.join(['C:\HE52_data\SIOPs\H2OabDefaults_FRESHwater.txt',
                           'dummyastar.txt',
                           'dummyastar.txt',
                           'C:\HE52_data\SIOPs\\astarphy_dinos_txt\\astarphy_dino_{0}.txt\n'.format(int(dinoD))])

    # LINE 16-17 - RECORD 5E: SPECIFIC SCATTERING PARAMTERS
    line16_17 = np.array([[0, -999, -999, -999, -999, -999],
                          [-666, -999, -999, -999, -999, -999]]).reshape(2, -1)

    # LINE 18-19 - RECORD 5F: SPECIFIC SCATTERING DATA FILE NAMES
    line18_19 = '\n'.join(['bstarDummy.txt',
                           'dummybstar.txt\n'])

    # LINE 20-21 - RECORD 5G: TYPES OF CONCENTRATIONS AND PHASE FUNCTIONS
    line20 = np.array([[0, 0, 550, 0.01, 0]]).reshape(1, -1)
    line21 = np.array([-2, 0, 0, 0, 0]).reshape(1, -1)

    # LINE 22-23 - RECORD 5H: PHASE FUNCITON FILE NAMES
    line22_23 = '\n'.join(['pureh2o.dpf',
                           'isotrop.dpf\n'])

    # LINE 24-75 - RECORD 6: WAVELENGTHS
    # 400:900 in 1 nm intervals = 500 'bands'
    line24 = np.array([85]).reshape(1, -1)
    line25_73 = np.arange(400, 800, 5).reshape(16, -1)
    line74 = np.arange(800,900,20).reshape(1,-1)
    line75 = np.array([900]).reshape(1, -1)

    # LINE 76 - RECORD 7: INELASTIC SCATTERING
    # AND LINE 77 - RECORD 8: SKY MODEL
    line76_77 = np.array([[0, 1, 0, 0, 4],
                          [2, 3, 30, 0, 0]]).reshape(2, -1)

    # REST OF LINES DO NOT NEED ATTENTION AND CAN STAY THE SAME
    # UNLESS YOU WANT TO CHANGE ATMOSPHERICS, DEPTHS, PATHS TO IOPS, ETC..
    line78 = np.array([-1, 0, 0, 29.92, 1, 80, 2.5, 15, 2, 300]).reshape(1, -1)  # atmospherics
    line79 = np.array([2, 1.34, 20, 35]).reshape(1, -1)  # surface info
    line80 = np.array([0, 0]).reshape(1, -1)  # bottom reflectance
    line81 = np.array([0, 5, 0, 1, 2, 3, 4, 5]).reshape(1, -1)  # depths

    # REST OF LINES ARE DATA FILE PATHS (82-93)...
    line82_93 = '\n'.join(['C:\HE52_data\SIOPs\H2OabDefaults_FRESHwater.txt',
                           '1',
                           'C:\HE52_data\NN\iopz_files\\acData' + name + '.txt',
                           'dummyFilteredAc9.txt',
                           'C:\HE52_data\NN\iopz_files\\bbData' + name + '.txt',
                           'C:\HE52_data\NN\chlz_files\Chlz' + name + '.txt',
                           'dummyCDOMdata.txt',
                           'dummyR.bot',
                           'dummydata.txt',
                           'dummyComp.txt',
                           'C:\HE52_data\NN\Ed_files\Ed_inputs.txt',
                           '..\data\MyBiolumData.txt'])

    #
    ''' --------- WRITE LINES TO IROOT FILE --------- '''

    with open(fname, 'w+') as f:
        np.savetxt(f, line1, fmt='%d,%d,%d,%1.3f,%d,%1.5f,%d,%1.1f')
        f.write(line2)
        f.write(line3)
        np.savetxt(f, line4, delimiter=',', fmt='%d')
        np.savetxt(f, line5, delimiter=',', fmt='%d')
        np.savetxt(f, line6, delimiter=',', fmt='%d')
        np.savetxt(f, line7, delimiter=',', fmt='%d')
        np.savetxt(f, line8_11, fmt='%d,%d,%d,%1.1f,%1.3f')
        f.writelines(line12_15)
        np.savetxt(f, line16_17, delimiter=',', fmt='%d')
        f.writelines(line18_19)
        np.savetxt(f, line20, fmt='%d,%d,%d,%1.2f,%d')
        np.savetxt(f, line21, delimiter=',', fmt='%d')
        f.writelines(line22_23)
        np.savetxt(f, line24, fmt='%d')
        np.savetxt(f, line25_73, delimiter=',', fmt='%d')
        np.savetxt(f, line74, delimiter=',',fmt='%d')
        np.savetxt(f, line75, fmt='%d')
        np.savetxt(f, line76_77, delimiter=',', fmt='%d')
        np.savetxt(f, line78, fmt='%d,%d,%d,%2.2f,%d,%d,%1.1f,%d,%d,%d')
        np.savetxt(f, line79, fmt='%d,%1.2f,%d,%d')
        np.savetxt(f, line80, delimiter=',', fmt='%d')
        np.savetxt(f, line81, delimiter=',', fmt='%d')
        f.writelines(line82_93)

        # WRITE FILE TO RUNLIST.TXT FOR HE5 BATCH PROCESSING
        with io.open(os.path.join(runlist, 'runlist.txt'), 'a+') as r:
            r.write(unicode(fname1)+'\r\n')

print '\n iRoot files generated!'



