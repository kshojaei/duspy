import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
E = np.linspace(0,30,2999)
Ei = np.linspace(0,30,26999)
path_dir = r'/Users/scienceman/Desktop/Project_Titan/Falcon/Dusty_Day_One/blackbox_photo3/'

T = pd.read_csv(path_dir + 'EEPF_4Pa_np5e15.txt')
TT = pd.read_csv(path_dir + 'EEPF_6Pa_np5e15.txt')
TTT = pd.read_csv(path_dir + 'EEPF_8Pa_np5e15.txt')
F = pd.read_csv(path_dir + 'EEPF_9Pa_np5e15.txt')
#FF = pd.read_csv(path_dir + 'EEPF_10Pa_dist_photo.txt')
#FFF = pd.read_csv(path_dir + 'EEPF_12Pa_dist_photo.txt')
#S = pd.read_csv(path_dir + 'EEPF_np5.5e14.txt')
EEPF_1 = np.array(T.iloc[:,1]).astype(float)
EEPF_2 = np.array(TT.iloc[:,1]).astype(float)
EEPF_3 = np.array(TTT.iloc[:,1]).astype(float)
EEPF_4 = np.array(F.iloc[:,1]).astype(float)
#EEPF_5 = np.array(FF.iloc[:,1]).astype(float)
#EEPF_6 = np.array(FFF.iloc[:,1]).astype(float)
#EEPF_7 = np.array(S.iloc[:,0]).astype(float)
#plt.title('Particle Density Variation')
fig = plt.figure()
#fig.set_size_inches(7,7)
import matplotlib.ticker as mticker

ax=plt.gca()
ax.plot(E, EEPF_1, linewidth=3, color= 'k')
ax.plot(E, EEPF_2,linewidth=3,color= 'b' )
ax.plot(E, EEPF_3,linewidth=3, color= 'r')
#plt.plot(E, EEPF_4,linewidth=3, color= 'r')
#plt.plot(E, EEPF_5, linewidth=3)
#plt.plot(E, EEPF_6, linewidth=3)
#plt.plot(E, EEPF_7, linewidth=4)


ax.grid(True)
#fig.set_size_inches(10, 10)

#plt.yticks(np.arange(0.015,0.03, step=0.0025))
#custom_ticks = np.linspace(0.01,0.03, 1, dtype=int)
#ax.set_yticks(custom_ticks)
ax.set_yscale('log')

from matplotlib.ticker import StrMethodFormatter, NullFormatter
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.02f}'))
ax.yaxis.set_minor_formatter(StrMethodFormatter('{x:0.02f}'))
#ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
plt.legend(['30 mTorr', '45 mTorr', '60 mTorr'])

ax.set_xlim(0,12)
plt.tick_params(axis="y",direction="in")
plt.rcParams.update({'font.size': 13})
plt.xlabel('Energy (eV)', fontsize=13)
plt.ylabel('$f_{o}$$(eV^{-3/2})$', fontsize=13)
plt.tick_params(direction='in')
plt.tick_params(which='minor', direction='in')
plt.rcParams.update({'font.size': 13})
ax.set_ylim(0.015,0.03)
#plt.text(0.3,0.02, r'$(a)$',fontsize=13)
#ax.set_yticks(np.arange(0.01,0.04, step=0.01))
#plt.text(0.5,0.037, r'(a)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#plt.legend(['30 mTorr','37.5 mTorr', '45 mTorr','60 mTorr', '75 mTorr', '90 mTorr', '105 mTorr'])
#plt.legend([])

plt.savefig('eepf_pow_paper.png', dpi = 600)


#print(np.trapz(EEPF_4 * np.sqrt(E), E))