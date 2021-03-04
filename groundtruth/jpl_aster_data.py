#%% 
'''
Organizes the JPL aster land emdmember library for use in 
adjacency RT modeling

Kravitz, 2021
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# path to data
path = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/ecospeclib-all_aster_jpl/'
#path = '/Users/jkravz311/Desktop/nasa_npp/groundtruth_data/test/'
#file = 'vegetation.tree.melaleuca.linariifolia.vswir.jpl198.jpl.asd.spectrum.txt'
def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))
files = listdir_nohidden(path)

# empty lists to append
infolist = []
spectralist = []

# spectra to ignore
bunk = ['mineral','tir','rock','meteorites','ancillary']


for file in files: 
    
    info = file.split('/')[-1]
    #print (info)
    info = info.split('.')
    
    # skip ignored spectra
    if any(x in info for x in bunk):   
        #print ('pass')
        pass
        
    
    else:  
        #print ('go')
        with open(file,'r') as f: 
            
            # get meta info for each spectra and store       
            line = f.readline()  
            name = line.split(': ')[1].rstrip() # Name
            line = f.readline()
            typ = line.split(': ')[1].rstrip() # Type
            line = f.readline() 
            clas = line.split(': ')[1].rstrip() # Class
            if typ in ['Vegetation','vegetation','non photosynthetic vegetation']:
                line = f.readline()
                genus = line.split(': ')[1].rstrip() # genus
                line = f.readline()
                species = line.split(': ')[1].rstrip() # species    
                subclass = np.nan # Subclass
                size = np.nan # Particle size
            else:
                line = f.readline()
                subclass = line.split(': ')[1].rstrip() # Subclass
                line = f.readline()
                size = line.split(': ')[1].rstrip() # Particle size   
                genus = np.nan # genus
                species = np.nan # species
            
            # continue skipping lines until spectra and stop
            while True:
                line = f.readline()
                if line in ['\n','\r\n']:
                    break
                
            # empty lists for lambdas and reflectances
            wl = []
            ref = []
            
            # read spectra data
            while True:
                line = f.readline()
                # stop at end of data
                if line in ['']:
                    break
                
                wl.append(line.split('\t')[0][1:].strip()) # some \t dont have space after
                ref.append(line.split('\t')[1].rstrip())                
                
            # create array of meta and append to list
            sampleData = np.array([name, typ, clas, genus, species, subclass, size])
            infolist.append(sampleData)
            
            # create spectra arrays and interp to 1nm resolution if not already
            newlambda = np.arange(.400,.900,.001)
            wl = np.array(wl).astype(float)
            ref = np.array(ref).astype(float)
            newref = np.interp(newlambda, wl, ref)
            spectralist.append(newref)
            
infodf = pd.DataFrame(infolist,columns=['Name','Type','Class','Genus','Species','Subclass','Particle size'])
spectradf = pd.DataFrame(spectralist, columns = np.arange(400,900))
finaldf = pd.concat([infodf,spectradf],axis=1)
finaldf.to_csv('/Users/jkravz311/Desktop/jpl_aster_data.csv',index=False)