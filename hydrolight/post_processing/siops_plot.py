##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

# plot params
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 18,
          'axes.titlesize': 16,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11
          }
pylab.rcParams.update(params)

# info
path = '/Users/jkravz311/Desktop/PhD/cyano_sims/hy_data/SIOPs/'

cyano = pd.read_csv(path+'cyano_siops.csv',index_col=0)
cyano = cyano.iloc[:-1,:]
dino_a = pd.read_csv(path+'dino_a.csv',index_col=0).T
dino_b = pd.read_csv(path+'dino_b.csv',index_col=0).T
dino_bb = pd.read_csv(path+'dino_bb.csv',index_col=0).T

wv = np.arange(400,905,5)

## PLOT -------------------------------------------------------

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3, figsize=(16,7))

# DINOS
# ABSORPTION
dino_a.iloc[:,4:45].plot(ax=ax1)
ax1.get_legend().remove()
#ax1.set_title('Specific Absorption',)
ax1.set_ylabel('a$^*$$_{phy}$ (m$^2mg^{-1})$',)
# SCATTERING
dino_b.iloc[:,4:45].plot(ax=ax2)
ax2.get_legend().remove()
#ax2.set_title('Specific Scatter',)
#ax2.set_xlabel('Wavelength (nm)',)
ax2.set_ylabel('b$^*$$_{phy}$ (m$^2mg^{-1})$',)
# BACKSCATTER
dino_bb.iloc[:,4:45].plot(ax=ax3)
ax3.get_legend().remove()
#ax3.set_title('Specific Backscatter',)
#ax3.set_xlabel('Wavelength (nm)',)
ax3.set_ylabel('bb$^*$$_{phy}$ (m$^2mg^{-1})$',)

# CYANOS
# ABSORPTION
cyano.a_aer.plot(ax=ax4, c='k', label='$\it{M. aeruginosa}$')
cyano.a_aph.plot(ax=ax4, c='r', label='$\it{Aphanizomenon\ flos-aquae}$')
cyano.a_ana.plot(ax=ax4, c='b', label='$\it{Anabaena\ cirinalis}$')
cyano.a_nod.plot(ax=ax4, c='g', label='$\it{Nodularia\ spumigena}$')
ax4.set_xlabel('Wavelength (nm)',)
ax4.set_ylabel('a$^*$$_{phy}$ (m$^2mg^{-1})$',)
ax4.legend(fontsize='small')
# scatter
cyano.b_aer.plot(ax=ax5, c='k', label='$\it{M. aeruginosa}$')
cyano.b_aph.plot(ax=ax5, c='r', label='$\it{Aphanizomenon\ flos-aquae}$')
cyano.b_ana.plot(ax=ax5, c='b', label='$\it{Anabaena\ cirinalis}$')
cyano.b_nod.plot(ax=ax5, c='g', label='$\it{Nodularia\ spumigena}$')
ax5.set_xlabel('Wavelength (nm)',)
ax5.set_ylabel('b$^*$$_{phy}$ (m$^2mg^{-1})$',)
# backscatter
cyano.bb_aer.plot(ax=ax6, c='k', label='$\it{M. aeruginosa}$')
cyano.bb_aph.plot(ax=ax6, c='r', label='$\it{Aphanizomenon\ flos-aquae}$')
cyano.bb_ana.plot(ax=ax6, c='b', label='$\it{Anabaena\ cirinalis}$')
cyano.bb_nod.plot(ax=ax6, c='g', label='$\it{Nodularia\ spumigena}$')
ax6.set_xlabel('Wavelength (nm)',)
ax6.set_ylabel('bb$^*$$_{phy}$ (m$^2mg^{-1})$',)
#
fig.tight_layout()
fig.savefig('/Users/jkravz311/Desktop/siop_plot.png',bbox_inches='tight',dpi=250)