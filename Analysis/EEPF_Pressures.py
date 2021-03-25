import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
E = np.linspace(0,30,2999)
Ei = np.linspace(0,30,26999)
path_dir = r'/Users/scienceman/Desktop/Project_Titan/Falcon/Dusty_Day_One/blackbox_photoemission2/'

T = pd.read_csv(path_dir + 'EEPF_4Pa_dist_photo.txt')
TT = pd.read_csv(path_dir + 'EEPF_5Pa_dist_photo.txt')
TTT = pd.read_csv(path_dir + 'EEPF_6Pa_dist_photo.txt')
F = pd.read_csv(path_dir + 'EEPF_8Pa_dist_photo.txt')
FF = pd.read_csv(path_dir + 'EEPF_10Pa_dist_photo.txt')
FFF = pd.read_csv(path_dir + 'EEPF_12Pa_dist_photo.txt')
#S = pd.read_csv(path_dir + 'EEPF_np5.5e14.txt')
EEPF_1 = np.array(T.iloc[:,1]).astype(float)
EEPF_2 = np.array(TT.iloc[:,1]).astype(float)
EEPF_3 = np.array(TTT.iloc[:,1]).astype(float)
EEPF_4 = np.array(F.iloc[:,1]).astype(float)
EEPF_5 = np.array(FF.iloc[:,1]).astype(float)
EEPF_6 = np.array(FFF.iloc[:,1]).astype(float)
#EEPF_7 = np.array(S.iloc[:,0]).astype(float)
#plt.title('Particle Density Variation')
fig = plt.figure()
import matplotlib.ticker as mticker

ax=plt.gca()
ax.plot(E, EEPF_1, linewidth=3, color= 'lime')
ax.plot(E, EEPF_2,linewidth=3,color= 'b' )
ax.plot(E, EEPF_3,linewidth=3, color= 'r')
#plt.plot(E, EEPF_4,linewidth=3, color= 'gold')
#plt.plot(E, EEPF_5, linewidth=3)
#plt.plot(E, EEPF_6, linewidth=3)
#plt.plot(E, EEPF_7, linewidth=4)


ax.grid(True)
#fig.set_size_inches(10, 10)
ax.set_ylim(0.01,0.045)

ax.set_yscale('log')

from matplotlib.ticker import StrMethodFormatter, NullFormatter
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.02f}'))
ax.yaxis.set_minor_formatter(StrMethodFormatter('{x:0.02f}'))
#ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())


ax.set_xlim(0,15)
plt.tick_params(axis="y",direction="in")
plt.rcParams.update({'font.size': 13})
plt.xlabel('Energy (eV)', fontsize=13)
plt.ylabel('$f_{o}$$(eV^{-3/2})$', fontsize=13)
plt.tick_params(direction='in')
plt.tick_params(which='minor', direction='in')
plt.rcParams.update({'font.size': 13})

plt.text(0.5,0.037, r'(a)', fontsize=13)
#plt.xticks(np.arange(0, 5, step=0.5))
plt.rcParams.update({'font.size': 13})
plt.legend(['30 mTorr','37.5 mTorr', '45 mTorr','60 mTorr', '75 mTorr', '90 mTorr', '105 mTorr'])


plt.savefig('eepf_pow_paper.png', dpi = 600)


#print(np.trapz(EEPF_4 * np.sqrt(E), E))
