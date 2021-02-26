from scipy import interpolate
import pandas as pd
import numpy as np

wl = np.arange(400,901,1)
idx440 = np.where(wl==440)

c_coef = pd.read_csv('/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/apc_C_coeficients_li_2015.csv',
                     index_col=0)
c1f = interpolate.interp1d(c_coef.index.values,c_coef.c1.values, fill_value='extrapolate')
c1new = c1f(wl)

c2f = interpolate.interp1d(c_coef.index.values,c_coef.c2.values, fill_value='extrapolate')
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
    astarpc = 0.008*admix ** -.451 # from Matthews 2013
    pc_li = apc620 / astarpc

    apc620_sim = aphy[220] - 0.24 * aphy[265]
    pc_sim = apc620_sim / astarpc

for k in range(6):
    pc_Li.append(round(pc_li,2))
    pc_simis.append(round(pc_sim,2))

fname = 'PC_' + nii+'.csv'

frames = [pc_Li,pc_simis]
pcData = np.vstack(frames).T

pcData = pd.DataFrame(pcData,index=depth[:,0],columns=['pc_Li','pc_simis'])
pcData.to_csv(os.path.join(pczpath,fname))