import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as co
from sympy.solvers import solve
from sympy import Symbol
import math
from scipy.optimize import fsolve
rp = 16E-9



ne = np.array([6.18E+14, 6.55E+14, 6.46E+14, 6.96E+14, 7.87E+14, 8.58E+14])
ni = np.array([1.33E+16, 1.60E+16, 1.82E+16, 2.17E+16, 2.37E+16, 2.43E+16])
floating = np.array([-1.14E+00, -1.38E+00, -1.58E+00, -1.89E+00, -2.06E+00, -2.11E+00])
charge = np.array([1.27E+01, 1.54E+01, 1.75E+01, 2.10E+01, 2.29E+01, 2.34E+01])
te = np.array([6.35E+00, 6.11E+00, 5.94E+00, 5.58E+00, 5.24E+00, 5.03E+00])
p = np.array([4,5,6,8,10,12])




ni = 2.43E+16
ne = 8.58E+14
#pressure = 12
Ti = 300
mass_one = 6.6335209E-26   
#mass_one = 3* 1.6735575E-27
Te = 5.03E+00
#pressure = np.array([2.5,4,6, 8])
pressure = np.array([30,37.5, 45,60,75,90])




def f(s):
    S = 4 * np.pi * (rp ** 2)
    A1 = S * ni * (((co.k*Ti)/(2* np.pi * mass_one))** 0.5)
    A2 = S * ne * (((co.k*Te*11604.52500617)/(2* np.pi * co.m_e))** 0.5)
    ve = A2 * np.exp((s*co.e)/(co.k*Te*11604.52500617))
    vi = A1 * (1-(s*co.e)/(co.k*Ti))
    return (ve - vi)

#print(fsolve(f,-1))

floating = np.array([-1.14E+00, -1.38E+00, -1.58E+00, -1.89E+00, -2.06E+00, -2.11E+00])
floating_oml = np.array([-3.09475712,-2.76608417,-2.45634849,-2.19311231,-2.15866484,-2.1950074])

print(len(floating), len(pressure))
#plt.plot(pressure, -floating)
#plt.plot(pressure, -floating_oml)


#plt.plot(pressure, (ne/ni))

fig = plt.gcf()
plt.grid(True)
plt.tick_params(direction= 'in')
plt.tick_params(axis="y",direction="in")
plt.text(31,2.1, r'(b)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#fig.set_size_inches(12, 12)
plt.plot(pressure, -floating, 'b', linewidth=4)
plt.plot(pressure, -floating_oml, 'b--', linewidth=4)

plt.ylabel('Particle Potential (eV)',fontsize=13)
plt.xlabel('Pressure (mTorr)', fontsize=13)
plt.legend(['With Secondary Emission Processes', 'OML Theory'], fontsize=13)
plt.ylim(1.0,3.3)
plt.xticks(np.arange(30,100, step=10))


#plt.ylim(1.6,3)
#plt.yticks(np.arange(1.6, 2.851, step=0.25))
plt.savefig('potential_paper_photo.png', dpi = 600)
#plt.savefig('potential_final.png', dpi = 600)

