import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as co
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator



#npp = np.array(T.iloc[0:8,0]).astype(float)

p = np.array([30,45,60,75,90,105,112.5])
temp = np.array([2.34E+03, 2.28E+03, 2.22E+03, 2.15E+03, 2.07E+03, 1.93E+03, 1.83E+03])
edea = np.array([8.63E-01, 7.05E-01, 5.25E-01, 3.38E-01, 1.64E-01, 5.33E-02, 3.22E-02])
edii = np.array([1.59E+02, 6.56E+01, 3.24E+01, 1.59E+01, 6.23E+00, 1.72E+00, 9.65E-01])



tot_heat = np.array([5.34E-11, 5.38E-11, 5.46E-11, 5.55E-11, 5.55E-11, 5.42E-11, 5.26E-11])

rad = np.array([4.23E+01, 3.30E+01, 2.54E+01, 1.92E+01, 1.34E+01, 7.74E+00, 5.14E+00])*(tot_heat/100)
cond = np.array([3.29E+01, 4.75E+01, 6.06E+01, 7.20E+01, 8.24E+01, 9.09E+01, 9.40E+01])*(tot_heat/100)
therm = np.array([2.47E+01, 1.95E+01, 1.40E+01, 8.75E+00, 4.15E+00, 1.35E+00, 8.26E-01])*(tot_heat/100)


eke = np.array([9.63E+01, 9.44E+01, 9.28E+01,9.15E+01, 9.03E+01, 8.90E+01, 8.83E+01])*(tot_heat/100)
ike = np.array([1.37E-01, 2.72E-01, 4.10E-01, 5.38E-01, 6.62E-01, 7.79E-01, 8.35E-01])*(tot_heat/100)
h3asdi = np.array([5.40E-01, 1.05E+00, 1.55E+00, 2.00E+00, 2.44E+00, 2.91E+00, 3.17E+00])*(tot_heat/100)
atomh = np.array([3.06E+00, 4.31E+00, 5.24E+00, 5.93E+00, 6.57E+00, 7.27E+00, 7.71E+00])*(tot_heat/100)
#ar3p2 = np.array([])*(tot_heat/100)
#ar3p0 = np.array([])*(tot_heat/100)




fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.grid(True)
plt.tick_params(direction= 'in')
#plt.tick_params(axis="y",direction="in")
#plt.text(31,2.4, r'(b)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#fig.set_size_inches(12, 12)
plt.plot(p, cond,'b--')
plt.plot(p, cond, 'o', color='b',markersize=7)

plt.plot(p, rad, 'k--')
plt.plot(p, rad, 'o', color= 'k',  markersize=7)


plt.plot(p, h3asdi,'r--')
plt.plot(p, h3asdi,'o', color= 'r',markersize=7)

plt.plot(p,atomh, '--', color= 'lime')
plt.plot(p,atomh, 'o', color= 'lime', markersize=7)


plt.plot(p, therm, '--', color= 'darkorange')
plt.plot(p, therm, 'o', color= 'darkorange', markersize=7)


plt.plot(p, eke, '--', color= 'fuchsia')
plt.plot(p, eke, 'o', color= 'fuchsia', markersize=7)
#plt.plot(p, eke, linewidth=3)


plt.plot(p, ike, '--', color= 'blueviolet')
plt.plot(p, ike,'o', color= 'blueviolet', markersize=7)

#plt.plot(temp,ar3p2, linewidth=3)
#plt.plot(temp, ar3p0, linewidth=3)

plt.xlim(27,92)
plt.xlabel('Pressure (mTorr)',fontsize=13)
plt.ylabel('Heat (W)', fontsize=13)
#plt.legend(['$Q_{rad}$','$Q_{Cond}$','$Q_{Therm}$','$Q_{EKE}$','$Q_{IKE}$','$Q_{H3+ Asso/Disso}$'
#            ,'$Q_{Atomic H}$'], fontsize=13)
#plt.xlim(15,71)
#plt.savefig('potential_paper.png', dpi = 600)
#plt.ylim(1.6,3)
#plt.yticks(np.arange(1.6, 2.851, step=0.25))
#plt.text(30,1.40e-11, r'$(d)$',color= 'k', fontsize=13)
plt.text(85,3.5e-11, r'$Q_{Cond}$',color= 'b', fontsize=13)
plt.text(85,0.65e-11, r'$Q_{rec}$',color= 'r',  fontsize=13)
plt.text(85,0.45e-11, r'$Q_{h}$',color= 'lime',  fontsize=13)
plt.text(85,0.16e-11, r'$Q_{ion}$',color= 'blueviolet',  fontsize=13)
plt.text(85,4.7e-11, r'$Q_{e}$',color= 'fuchsia',  fontsize=13)
plt.text(85,1.5e-11, r'$Q_{rad}$',color= 'k',  fontsize=13)
plt.text(64,0.8e-11, r'$Q_{emission}$',color= 'darkorange',  fontsize=13)


'''
plt.text(46, 6.1, r'(c)', fontsize=13)
plt.text(80, 5, r'$T_{e}$', color= 'k', fontweight= 'bold', fontsize=13)
#plt.yticks(np.arange(5, 7.5, step=1))
plt.tick_params(which='minor',color='k')
plt.minorticks_on()
plt.tick_params(axis='x', which='minor', bottom=False)
plt.xlim(45,100)
plt.ylim(4.5,6.5)

plt.savefig('plot_temperature.png', dpi = 600)

'''



