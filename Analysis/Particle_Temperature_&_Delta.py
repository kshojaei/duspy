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
thii = np.array([6.61E+00, 3.44E+00, 1.85E+00, 4.12E-01, 3.07E-02, 6.95E-04])
photoii = np.array([6.83E-01, 4.78E-01, 3.70E-01, 2.74E-01, 2.38E-01, 2.28E-01])
quenii = np.array([1.08E-02, 9.09E-03, 8.16E-03, 7.54E-03, 7.64E-03, 8.26E-03])


fig, ax1 = plt.subplots()


#ax1.rcParams.update({'font.size' : 13})

color = 'tab:red'
ax1.set_xlabel('Pressure (mTorr)', fontsize=13)
ax1.set_ylabel('${T_p} (K)$', color=color, fontsize=13)
ax1.tick_params(which='minor', direction='in')
ax1.tick_params(direction= 'in')
ax1.plot(p,temp, color='r', linewidth=3, marker='o')
ax1.set_ylim(990,1820)
plt.xticks(np.arange(30,120, step=10))
plt.xlim(28,92)
#ax1.tick_params(axis='y', labelcolor=color)
plt.grid(True)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax1.text(31,1220, r'$(c)$',fontsize=13)
color = 'tab:blue'
ax2.set_ylabel('$\delta_{T}$', color='fuchsia', fontsize=13)  # we already handled the x-label with ax1
ax2.tick_params(which='minor', direction='in')
ax2.tick_params(direction= 'in')
ax2.plot(p,edii, color='fuchsia', linewidth=3, marker= 'o')
#ax2.plot(p,thii, color='lime', linewidth=1, marker= 's')
#ax2.plot(p,photoii, color='g', linewidth=1, marker= '^')
#ax2.plot(p,qii, color='b', linewidth=3, marker= 'o')
ax2.set_ylim(-1,8.5)
#ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
fig.savefig('temp_prop.png', dpi = 600)


'''
fig = plt.gcf()
plt.grid(True)
plt.tick_params(direction= 'in')
plt.tick_params(axis="y",direction="in")
#plt.text(31,2.4, r'(b)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#fig.set_size_inches(12, 12)
#plt.plot(p, edii, 'b', linewidth=3, marker= 'o')
#plt.plot(p,thii, 'r', linewidth=3, marker= 's')
#plt.plot(p,photoii, 'lime', linewidth=3, marker= '^')
plt.plot(p, thii/edii*100, color= 'lawngreen', linewidth='3', marker= 's')
plt.plot(p, photoii/edii*100, color= 'fuchsia', linewidth='3', marker= 'o')
plt.ylabel('Percentage (%)')
plt.legend(['$\delta_{therm}$ / $\delta_{T}$','$\delta_{photo}$ / $\delta_{T}$' ])
#plt.plot(p,quenii, 'b--', linewidth=4)

#plt.ylabel('Particle Potential (eV)',fontsize=13)
plt.xlabel('Pressure (mTorr)', fontsize=13)
#plt.legend(['ed / ii', 'therm / ii', 'photo / ii'], fontsize=13)
#plt.ylim(0,8)
plt.xticks(np.arange(30,94, step=10))

'''
'''
ne = np.array([6.18E+14, 6.55E+14, 6.46E+14, 6.96E+14, 7.87E+14, 8.58E+14])
ni = np.array([1.33E+16, 1.60E+16, 1.82E+16, 2.17E+16, 2.37E+16, 2.43E+16])
floating = np.array([-1.14E+00, -1.38E+00, -1.58E+00, -1.89E+00, -2.06E+00, -2.11E+00])
charge = np.array([1.27E+01, 1.54E+01, 1.75E+01, 2.10E+01, 2.29E+01, 2.34E+01])
te = np.array([6.35E+00, 6.11E+00, 5.94E+00, 5.58E+00, 5.24E+00, 5.03E+00])
#p = np.array([4,5,6,8,10,12])


#p = np.array([4,6,8,10,12, 14])

fig = plt.gcf()
plt.grid(True)
plt.tick_params(direction= 'in')
plt.tick_params(axis="y",direction="in")
#plt.text(31,2.4, r'(b)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#fig.set_size_inches(12, 12)
#plt.plot(p, te, 'r', linewidth=3)
plt.plot(p,ne*1e-6, 'b', linewidth=3)
plt.plot(p,ni*1e-6, 'r', linewidth=3)
#plt.plot(p,quenii, 'b--', linewidth=4)
plt.yscale('log')
#plt.ylim(1e8, 1e11)
#
plt.ylabel(' Density (1/cm3)',fontsize=13)
plt.xlabel('Pressure (mTorr)', fontsize=13)
plt.legend(['ne', 'ni'])
#plt.legend(['ed / ii', 'therm / ii', 'photo / ii'], fontsize=13)
#plt.xlim(28,112)
plt.xticks(np.arange(30,95, step=10))
'''







