import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as co
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator



#npp = np.array(T.iloc[0:8,0]).astype(float)

p = np.array([30,37.5, 45,60,75,90])
temp = np.array([1.79E+03, 1.66E+03, 1.55E+03, 1.35E+03, 1.17E+03, 1.02E+03])
edea = np.array([7.45E-01, 6.24E-01, 4.82E-01, 2.37E-01, 1.23E-01, 1.16E-01])
edii = np.array([7.34E+00, 3.93E+00, 2.23E+00, 6.91E-01, 2.75E-01, 2.36E-01])

tot_heat = np.array([1.82E-11, 1.82E-11, 1.85E-11, 1.91E-11, 1.92E-11, 1.89E-11])

rad = np.array([1.27E+01, 7.37E+00, 4.53E+00, 1.80E+00, 7.22E-01,3.01E-01])*(tot_heat/100)
cond = np.array([7.08E+01, 8.07E+01, 8.74E+01, 9.51E+01, 9.79E+01, 9.85E+01])*(tot_heat/100)
therm = np.array([1.64E+01, 1.19E+01, 8.10E+00, 3.13E+00, 1.39E+00, 1.25E+00])*(tot_heat/100)


eke = np.array([7.69E+01, 6.90E+01, 6.25E+01, 5.12E+01, 4.38E+01, 4.07E+01])*(tot_heat/100)
ike = np.array([2.65E+00, 3.98E+00, 5.17E+00, 7.35E+00, 8.71E+00, 9.15E+00])*(tot_heat/100)
h3asdi = np.array([9.76E+00, 1.42E+01, 1.80E+01, 2.49E+01, 2.95E+01, 3.14E+01])*(tot_heat/100)
atomh = np.array([1.07E+01, 1.28E+01, 1.43E+01, 1.66E+01, 1.80E+01, 1.87E+01])*(tot_heat/100)
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
plt.text(85,1.7e-11, r'$Q_{Cond}$',color= 'b', fontsize=13)
plt.text(85,0.7e-11, r'$Q_{rec}$',color= 'r',  fontsize=13)
plt.text(85,0.4e-11, r'$Q_{h}$',color= 'lime',  fontsize=13)
plt.text(85,0.21e-11, r'$Q_{ion}$',color= 'blueviolet',  fontsize=13)
plt.text(85,0.85e-11, r'$Q_{e}$',color= 'fuchsia',  fontsize=13)
plt.text(75,0.21e-11, r'$Q_{rad}$',color= 'k',  fontsize=13)
plt.text(65,0.21e-11, r'$Q_{emission}$',color= 'darkorange',  fontsize=13)


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



