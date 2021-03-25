import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as co
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator

#path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/01-31-Cu_double_Ar/Ar Addition/ne_ni_Te_hu.txt'
path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/12-19-copper/pressure study/ne_ni_Te_hu.txt'
T = pd.read_csv(path_dir, delimiter= ',')

'''
old ar as the dominant ion wrong
#npp = np.array(T.iloc[0:8,0]).astype(float)
ni = np.array([5.4182850127153837e+17,5.575845067809493e+17,5.644985321522464e+17,
               5.633676962890809e+17,5.59653141903224e+17,5.579332919403409e+17])*1e-6
ne = np.array([1.279906940691818e+16,1.3606553138772766e+16,1.3886957663548376e+16,
               1.408202574044621e+16,1.4046611190245878e+16,1.3949349097469878e+16])*1e-6*4.386389364536861
Te = np.array([4.97463708631586,5.098855971950952,5.273347310076409,
               5.470563068244049,5.706148238154203,5.966881334713291])
Pin = np.array([46.5, 56.7, 66.1, 75.5, 84.6, 95.5])

'''

#npp = np.array(T.iloc[0:8,0]).astype(float)
ni = np.array([1.4908605625747762e+16,1.533979633219836e+16,1.5530744792708636e+16,
               1.5501790202149272e+16,1.5398119011684332e+16,1.533979633219836e+16])*1e-6


ne = np.array([1281092636190372.8,1360655313877276.8, 1389059203807093.5,
               1409629227797907.2,1405337688753337.5,1360655313877276.8])*1e-6


Te = np.array([4.977605318945703,5.09885597200194,5.274358098095077,5.474586082601373,
               5.708184836696361,5.972755229466758])
Pin = np.array([46.5, 56.7, 66.1, 75.5, 84.6, 95.5])



Te_unc = np.array([0.11806379246428325,0.10519404202760535,0.14505369970578208,
                   0.12072011376614825,0.1243792967107556,0.11939424989503714])*Te

ne_unc = np.array([0.1373651757339936,0.13374970848488985,0.19003259554513072,
                   0.17259288635459133,0.17840062614987706,0.16953358222440335])*ne

ni_unc = np.array([0.06998951839618957,0.06946813230944406,0.06819520813671882,
                   0.06747619183155132,0.06693515754603287,0.06945630664656195])*ni




#print(ni)


print(ni/ne)
plt.subplot(2,1,1)
plt.plot(Pin, ni,'r--',linewidth=3)
plt.errorbar(Pin, ni, xerr=0.2, yerr=ni_unc, fmt='ro')
plt.plot(Pin, ne,'b--', linewidth=3)
plt.errorbar(Pin, ne, xerr=0.2, yerr=ne_unc, fmt='bo')
#plt.legend(['$n_{i}$','$n_{e}$'])
plt.tick_params(direction= 'in')
plt.xticks(Pin, " ")
plt.rcParams.update({'font.size' : 13})
plt.tick_params(axis = "x", which = "both", bottom = False, top = False)
plt.ylabel('n ($cm^{-3})$', fontsize=13)
plt.yscale('log')
plt.tick_params(which='minor', direction='in')
plt.xlim(45,100)
plt.text(46, 5e9, r'(b)', fontsize=13)
plt.text(80, 0.8e10, '$n_{i}$',color= 'r', fontsize=13, fontweight= 'bold')
plt.text(80, 2e9, '$n_{e}$',color= 'b', fontsize=13, fontweight= 'bold')
plt.ylim(1e9,3e10)
plt.subplot(2,1,2)
plt.plot(Pin, Te, 'k--',linewidth=3)
plt.tick_params(which='minor', direction='in')
plt.errorbar(Pin, Te, xerr=0.2, yerr=Te_unc, fmt='ko')
plt.tick_params(direction= 'in')
plt.ylabel('$T_{e}$(eV)', fontsize=13)
plt.xlabel('Pressure (mTorr)', fontsize=13)
plt.rcParams.update({'font.size' : 13})
plt.text(46, 6.1, r'(c)', fontsize=13)
plt.text(80, 5, r'$T_{e}$', color= 'k', fontweight= 'bold', fontsize=13)
#plt.yticks(np.arange(5, 7.5, step=1))
plt.tick_params(which='minor',color='k')
plt.minorticks_on()
plt.tick_params(axis='x', which='minor', bottom=False)
plt.xlim(45,100)
plt.ylim(4.1,7)

plt.savefig('plot_temperature.png', dpi = 600)





