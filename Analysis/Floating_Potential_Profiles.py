import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as co
from sympy.solvers import solve
from sympy import Symbol
import math
from scipy.optimize import fsolve
rp = 16E-9



ne = np.array([1.53E+14, 1.41E+14, 1.32E+14, 1.25E+14, 1.22E+14, 1.27E+14, 1.33E+14])
ni = np.array([1.12E+16, 1.59E+16, 1.96E+16, 2.25E+16, 2.50E+16, 2.70E+16, 2.78E+16])
floating = np.array([-1.99E-01, -2.84E-01, -3.50E-01,-4.03E-01, -4.48E-01, -4.84E-01,-4.98E-01])
charge = np.array([2.21E+00, 3.15E+00, 3.89E+00, 4.48E+00, 4.97E+00, 5.37E+00, 5.53E+00])
te = np.array([7.68E+00, 7.67E+00, 7.69E+00, 7.74E+00, 7.70E+00, 7.54E+00,7.35E+00])
p = np.array([4,6,8,10,12,14,15])




ni = 2.78E+16
ne = 1.33E+14
#pressure = 12
Ti = 300
mass_one = 6.6335209E-26   
#mass_one = 3* 1.6735575E-27
Te = 7.35E+00
#pressure = np.array([2.5,4,6, 8])
pressure = np.array([30,45,60,75,90,105,112.5])




def f(s):
    S = 4 * np.pi * (rp ** 2)
    A1 = S * ni * (((co.k*Ti)/(2* np.pi * mass_one))** 0.5)
    A2 = S * ne * (((co.k*Te*11604.52500617)/(2* np.pi * co.m_e))** 0.5)
    ve = A2 * np.exp((s*co.e)/(co.k*Te*11604.52500617))
    vi = A1 * (1-(s*co.e)/(co.k*Ti))
    return (ve - vi)

print(fsolve(f,-1))

floating = np.array([-1.99E-01, -2.84E-01, -3.50E-01,-4.03E-01, -4.48E-01, -4.84E-01,-4.98E-01])
floating_oml = np.array([-1.35165658,-0.91937701,-0.71273854,-0.595131,-0.52310208,-0.49868573,-0.49990725])


print(len(floating), len(pressure))
#plt.plot(pressure, -floating)
#plt.plot(pressure, -floating_oml)


#plt.plot(pressure, (ne/ni))

fig = plt.gcf()
plt.grid(True)
plt.tick_params(direction= 'in')
plt.tick_params(axis="y",direction="in")
#plt.text(31,0.7, r'(b)', fontsize=13)
plt.rcParams.update({'font.size': 13})
#fig.set_size_inches(12, 12)
plt.plot(pressure, -floating, 'b', linewidth=4)
plt.plot(pressure, -floating_oml, 'b--', linewidth=4)

plt.ylabel('Particle Potential (eV)',fontsize=13)
plt.xlabel('Pressure (mTorr)', fontsize=13)
plt.legend(['With Secondary Emission Processes', 'OML Theory'], fontsize=13)
#plt.ylim(1.0,3.3)
plt.xticks(np.arange(30,120, step=10))
plt.xlim(28,115)

#plt.ylim(1.6,3)
#plt.yticks(np.arange(1.6, 2.851, step=0.25))
plt.savefig('potential_paper_photo.png', dpi = 600)
#plt.savefig('potential_final.png', dpi = 600)

