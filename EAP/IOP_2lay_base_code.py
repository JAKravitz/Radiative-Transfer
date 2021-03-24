#%%
# fortran wrapper
from numpy import f2py
sourcefile = open('/Users/jkravz311/git_projects/Radiative-Transfer/EAP/Dmmex_R14B_4.f','rb')
sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='Dmmex_R14B_4')

import Dmmex_R14B_4
import scipy.io as io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
l = np.arange(.4, .905, .005) # wavelength range and resolution (changing this changes your interp value when normalising kshell)
int_val = 55 # refers to l[55] which is 675 nm. 255 for 1 nm resolution

Vs = 0.2
ci = 2.5e6
D_eff = np.array([15])

V_eff= 0.6
mf = io.loadmat('/Users/jkravz311/git_projects/Radiative-Transfer/EAP/501nm_extended_e1701000.mat')

# using absorption for imaginary refractive index
im = pd.read_csv('/Users/jkravz311/Desktop/phyto_imag_refr_idx.csv')
im = im.filter(regex='^[0-9]')
im_wv = np.arange(.4,.801,.001)
im_a1 = im.iloc[0,:].values


Vc=1 - Vs
FR=(1- Vs) ** (1/ 3)# relative volume of core to shell
nmedia = 1.334
wavelength = l/ nmedia

wvno = 2 * np.pi / wavelength # this is a Mie param - combo of size and wavelength

# hilbert transform
def analytic_signal(x):
    from scipy.fftpack import fft, ifft
    N = len(x)
    X = fft(x, N)
    h = np.zeros(N)
    h[0] = 1
    h[1:N//2] = 2* np.ones(N// 2-1)
    h[N// 2] = 1
    Z = X * h
    z = ifft(Z, N)
    return z

from scipy.interpolate import griddata
kcore = griddata(mf['RIs'][:, 5], mf['RIs'][:, 0], l, 'linear')
#kshell_base = griddata(mf['RIs'][:, 5], mf['RIs'][:, 2], l, 'linear')

# measured abs
kshell_base = griddata(im_wv, im_a1, l, 'linear',)
# replace negative/zeros with nans
kshell_base = np.where(kshell_base < 0, np.nan, kshell_base)
kshell_base = np.where(kshell_base == 0, np.nan, kshell_base)
# replace nans with minimum value
kshell_base[np.isnan(kshell_base)] = min(kshell_base)

kshell_norm = (6.75e-7/ nmedia) * (0.027 * ci/ Vs) / (4 * np.pi) #scale to this theoretical max unpackaged chl abs at 675 nm
kshell = kshell_base * (kshell_norm / kshell_base[int_val]) #55 is index of 675 nm

nshell = 1.14 + np.imag(analytic_signal(kshell))
ncore = 1.02 + np.imag(analytic_signal(kcore))
khom = kcore*Vc + kshell*Vs # real refractive index 
nhom = ncore*Vc + nshell*Vs
mshell = nshell - kshell*1j
mcore = ncore - kcore*1j
mhom = nhom - khom*1j

psd = np.arange(1,101,1)
deltad=1
deltadm=1
theta901 = np.arange(0, 90.1, 0.1) # length 901

nang=901
angles=np.cos(np.deg2rad(theta901)) # length 901
theta1 = np.deg2rad(theta901) 
theta2 = np.flipud(np.pi-np.deg2rad(theta901[0:900]))
theta=np.append(theta1,theta2)

back1=np.where(theta==(np.pi/2)) 
back2=np.where(theta==(np.pi))
d1=np.diff(theta)
dtheta = np.append(d1, d1[0])

# preparing variables to be filled
VSF = np.zeros((len(D_eff),len(l), 1801))   #dimensions jjj(deff), nii(wavelength), 1801 angles    
PF_check = np.zeros((len(D_eff),len(l)))
d_alpha = []  
PF = np.zeros((len(D_eff), len(wavelength), 1801))

#%%
# declare all lists and arrays to fill in the jjj loop (refilled each iteration)
Qc, Sigma_c, c, Qb, Sigma_b, b, Qa, Sigma_a, a, Qbb, Sigma_bb, bb, bbtilde = (np.zeros([len(D_eff),len(l)]) for i in range(13))
a_sol, a_solm, Qastar2_dir, Qas21 = (np.zeros([len(D_eff),len(l)]) for i in range(4))

for nii in np.arange(0,len(l)): # this is the wavelength loop
    print(nii)

    # declare lists to be filled on each iteration of the psd loop
    II, phaseMB, alpha, bbprob, bbprob1, Qbbro, checkMB, Qbro, Qcro, M1 = ([] for i in range (10))
    

    for jj in np.arange(0, len(psd)): # this is the psd loop

        [QEXT,QSCA,GQSC,QBS,m1,m2,s21,d21] = Dmmex_R14B_4.dmilay((psd[jj]*FR)/2,psd[jj]/2,wvno[nii],mshell[nii],mcore[nii],angles,901,901)
    
        # on each iteration of jj, we get a different QEXT and QSCA out. So these must be stored in their own array
        Qcro.insert(jj,QEXT) 
        Qbro.insert(jj,QSCA)    
        
        m1_seta = [num[0] for num in m1]
        m1_setb = [num[1] for num in m1]
        M1 = np.append(m1_seta,m1_setb[900:0:-1])
        M2 = np.append(m2[0:901,0], m2[900:0:-1,1])
        myval = (M1+M2)/2 
        II.insert(jj, myval) 
        
        alpha2=2*np.pi*(psd[jj]/2)/wavelength[nii] 
        alpha.insert(jj, alpha2) 

        phaseMB_jj = [II[jj] / (np.pi* Qbro[jj]* (alpha[jj]**2))]
        phaseMB.insert(jj,phaseMB_jj)

        checkMB_jj = [2* np.pi* np.sum(phaseMB_jj * np.sin(theta) * dtheta)]
        checkMB.insert(jj,checkMB_jj)

        section_jj = [item[900:1801] for item in phaseMB_jj]
        bbprob_jj = 2*np.pi* np.sum((section_jj *np.sin(theta[900:1801]) *dtheta[900:1801]))
        bbprob.insert(jj, bbprob_jj) 
        Qbbro_jj = QSCA * bbprob_jj 
        Qbbro.insert(jj,Qbbro_jj) 
    
	# we are still in the nii loop here! just the jj loop has ended

    d_alpha_nii = alpha[1] - alpha[0]
    d_alpha.insert(nii,d_alpha_nii)

	# jjj loop starts here
    for jjj in np.arange(0,len(D_eff)):

        exponent = (-psd/ 2)/ ((D_eff[jjj]/ 2) * V_eff)
        psd2 = 1.0e20 * np.power((psd/2),((1-3* V_eff)/V_eff)) * np.exp(exponent)
        psdm1 = psd / 1e6; psdm2 = psd2 * 1e3; civol = np.pi/ 6 * sum(psdm2 * psdm1 **3 * deltadm)
        psdm2 = psdm2 * (1./ (civol * ci))
        psdvol = np.pi/6 * sum(psdm2 * np.power(psdm1, 3) * deltadm)
		
        Qc[jjj, nii] = sum(Qcro *psdm2 * np.power(psdm1,2) * deltadm)/ sum(psdm2 * np.power(psdm1,2) *deltadm)
        Sigma_c[jjj,nii] = np.pi/4 * Qc[jjj, nii] * sum(np.power(psdm1, 2) * deltadm)
        c[jjj,nii] = np.pi/4* Qc[jjj, nii]* sum(psdm2* np.power(psdm1,2)* deltadm)
        
        Qb[jjj, nii] = sum(Qbro * psdm2 * np.power(psdm1,2) * deltadm) /sum(psdm2* np.power(psdm1,2)* deltadm) 	            
        Sigma_b[jjj,nii] = np.pi/4 * Qb[jjj,nii]* sum(np.power(psdm1,2)* deltadm)
        b[jjj, nii] = np.pi/4* Qb[jjj, nii]* sum(psdm2* np.power(psdm1,2)* deltadm)

        Qbb[jjj, nii] = sum(Qbbro * psdm2 * np.power(psdm1,2) * deltadm) /sum(psdm2* np.power(psdm1,2)* deltadm)
        Sigma_bb[jjj, nii] = np.pi/4 * Qbb[jjj, nii] * sum(np.power(psdm1, 2) * deltadm)
        bb[jjj, nii] =  np.pi/4* Qbb[jjj, nii]* sum(psdm2 * np.power(psdm1, 2) * deltadm)
        
        Qa[jjj, nii] = Qc[jjj, nii] - Qb[jjj, nii]
        Sigma_a[jjj, nii] = np.pi/4 * Qa[jjj, nii]* sum(np.power(psdm1,2)* deltadm)
        a[jjj, nii] = c[jjj, nii] - b[jjj, nii]

        betabar, VSF_1 = ([] for i in range(2))
        checkbar = []
        b_check, bb_check = (np.zeros((len(D_eff),len(wavelength))) for i in range(2))
        
        bbtilde[jjj, nii] = bb[jjj, nii] / b[jjj, nii]
			
		# this little sub loop is INSIDE the jjj loop		
        for ai in range (0, nang * 2 - 1): # this should go 1801 times - doesn't include the last in range  
    		   # need a variable to get(II(:,ai)):
            varII = [item[ai] for item in II]
            betabar_ai = (1 / np.pi) * (sum(varII * psdm2 * d_alpha[nii]) / sum(Qbro * psdm2 * np.power(alpha, 2) * d_alpha[nii]))
            betabar.insert(ai, betabar_ai)
            VSF_1_ai = betabar[ai] * b[jjj, nii]  
            VSF_1.insert(ai, VSF_1_ai) # this gives VSF_1 of 1801 angles. For the current instance of nii (wavelength) and jjj (Deff)

        # checkbar is back outside of the sub loop
        checkbar = (2* np.pi * sum(betabar * np.sin(theta) * dtheta))
        PF_check[jjj,nii] = checkbar
        
        b_check[jjj,nii] = 2 * np.pi * sum((VSF_1) * np.sin(theta) * dtheta)
        
        PF[jjj,nii,:] = betabar
        VSF[jjj,nii,:] = VSF_1 # VSF_1s are put into matrix on each iteration of Deff, then wavelength.
        # We want to get out all wavelengths but backscatter only angles for each Deff:
            
        slice_obj = slice(900, 1801) 
        VSF_b = VSF[jjj, nii, slice_obj] # want to get the backward angles for this instance of Deff and all the wavelengths 
        bb_check[jjj,nii] = 2 * np.pi * sum((VSF_b) * np.sin(theta[900: 1801]) * dtheta[900: 1801])      
        
        ## Package effect calculations:
	
        ##Set up all phase lag and optical thickness parameters
        ##Calculate acm (absorption celular material) for all layers
        acm_core = (4 * np.pi * kcore[nii]) / (wavelength[nii] * 1e-6)
        acm_shell=(4 * np.pi * kshell[nii]) / (wavelength[nii] * 1e-6)
        acm_hom = (4 * np.pi * khom[nii]) / (wavelength[nii] * 1e-6)
        q = (D_eff[jjj] / 2 * FR) / (D_eff[jjj] / 2)

        Qas21[jjj, nii] = 3 / 2 * (Qa[jjj, nii] / (acm_core * np.power(q, 3) + acm_shell * (1 - np.power(q, 3)) * 1e-6 * D_eff[jjj]))
      
        ##Direct volume equivalent determination of package effect
        a_sol[jjj, nii] = psdvol * Vc * acm_core + psdvol * Vs * acm_shell
        a_solm[jjj, nii] = psdvol * acm_hom
        Qastar2_dir[jjj,nii] = a[jjj, nii] / a_sol[jjj, nii] #Hayley for input into fluorescence algorithm
      
 	    # both the jjj loop and the nii loop end here.
         
#%% plot results

a = pd.DataFrame(a,columns=l)
a.to_csv('/Users/jkravz311/Desktop/phyto_test/a_V0.4_C6_n1.14.csv')
b = pd.DataFrame(b,columns=l)
b.to_csv('/Users/jkravz311/Desktop/phyto_test/b_V0.4_C6_n1.14.csv')
bb = pd.DataFrame(bb,columns=l)
bb.to_csv('/Users/jkravz311/Desktop/phyto_test/bb_V0.4_C6_n1.14.csv')

a.T.plot()
b.T.plot()
bb.T.plot()




