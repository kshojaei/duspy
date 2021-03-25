from __future__ import division
import numpy as np 
import pandas as pd
from bolos import solver, grid, parser, process, target
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as co
import scipy as sp
import miepython
import matplotlib.ticker as mticker
'''############# Input Values ############## '''
ne_guess = 3.01e15                # Electron Density [1/m3] (Initial Guess)
Vp_guess = -0.9355485050377268# Particle Potential (eV) (Initial Guess)
ni_guess = 1.1878518915640684e+16 # Ion Density [1/m3]
EN_range = np.linspace(0,1,int(1E1))# Range of iterable Reduced Electric Field 
ni_range = np.linspace(0, -1E17,int(1E3))
amp_range = np.linspace(0,-1E-17,int(1E3))
rp = 16E-9                          # Particle Radius [m]
RR = 0.1                            # Maximum Discharge Dimension [m]
npp = 1E15                          # Particle Density [1/m3]
Te_guess = 5.6                      # Electron Temperature [eV]
EN_guess = 760            # Reduced Electric Field [Td]
p = 4                           # Pressure [Pa]
tol_iter = 0.01                     # Tol with Vp and EN
P_in = 90                           # Power Transferred into the Plasma [W]
f_maxmin = 0.9                      # Ratio of gases (Dominant/Subservient)
mass_one = 6.6335209E-26            # Mass of the Dominant Gas (Ar) [kg]
mass_two = 2* 1.6735575E-27         # Mass of the Subservient Gas H2 [kg]
mass_dom = 3* 1.6735575E-27         # Mass of the dominant ion H3+ [kg]
mass_H = 1.6735575E-27              # Hyrogen Atomic Mass [kg]
M_Ar = 39.948E-3                    # Atomic Weight of Argon  [kg/mol]
W_bulk = 4.7                        # Bulk Work Function [eV]
E_fermi = 7                         # Copper Fermi Energy (eV)
laser_pow = 0                       # Laser Power [W] Gaussian Profile: (Beam Diameter)
laser_spot_radius = 695E-6/2        # Laser Beam Radius [m]
laser_wavel = 532E-9                # Laser Wavelength [m]
mtwo = 1.06                         # beam quality factor 
Vsh = 1300                           # Sheath Voltage [V]
Q = 1                               # Quantum Yield 0.25-1.0 electron/ incident photon
R = 9E-2                            # Cylindrically Shaped Plasma Radius [m]
height = 10.16E-2                   # Cylindrically Shaped Plasma Height [m]
gold_dir = r'gold.txt'              # Data: Complex Index of Refraction vs. Wavelength [m]
copper_dir = r'copper.txt'          # Data: Complex Index of Refraction vs. Wavelength [m]
E = np.linspace(0,20,2000)          # Electron Energy Range [eV]
bb_guess = 1.42E-13               # Amplitude of detachment cross section
gs = grid.LinearGrid(0,30,3000)
T_gas = 300                         # Gas Temperature
Ti = 300                            # Ion Temperature [K]
tolerance = 1E-5                    # Convergence Tolerance 
# Particle Mass [kg]
m_particle = 8960 * ((4/3) * np.pi * (rp**3))
CS_dir = 'Ar_H2.txt'
gama_1 = 1.32E8                     # Transition Probability (3P1) [1/s]
gama_2 = 5.32E8                     # Transition Probability (1P1) [1/s]
f_1 = 0.0675                        # Oscillator Strength for 106.7 nm (3P1) Ex Argon
f_2 = 0.2629                        # Oscillator Strength for 104.8-9 nm (1P1) Ex Argon
# Escape factor for 104.8-9 nm (1P1) ex argon and 106.7 nm (3P1) ex argon 
# as a function of pressuree (mTorr)
# # Escape Factor (3P1) Ex Argon lmda_1: 106.7E-9 nm
def escape_factor_1p1(pressure):
    slope = (0.8E-3- 2.5E-3)/(100 - 30)
    rel = 0.8E-3 + slope*(pressure - 100)
    return rel
# Escape Factor (1P1) Ex Argon lmda_2: 104.8E-9 nm
def escape_factor_3p1(pressure):
    slope = (0.2E-3- 2.0E-3)/(100 - 30)
    rel = 0.2E-3 + slope*(pressure - 100)
    return rel
#g_escape_1 = escape_factor_1p1(30)        # Escape Factor (3P1) Ex Argon lmda_1: 106.7E-9 nm
#g_escape_2 = escape_factor_3p1(30)               # Escape Factor (1P1) Ex Argon lmda_2: 104.8E-9 nm
g_escape_1 = 2E-3               # Escape Factor (3P1) Ex Argon lmda_1: 106.7E-9 nm
g_escape_2 = 2.5E-3               # Escape Factor (1P1) Ex Argon lmda_2: 104.8E-9 nm
charge_lim = (4*np.pi*co.epsilon_0* W_bulk/co.e)*rp + 3/8
print('Charge Limit:' , charge_lim)
'''######################## Gas Cross Sections ########################### '''
# Gas Effective Cross Section [m2]
def Gas_Effective(E,species_name):
    with open(CS_dir) as fp:
        processes = parser.parse(fp)
        for i in processes:
            while (i['target']) == species_name and i['kind'] == 'EFFECTIVE':
                res = [list(z) for z in zip(*i['data'])]
                eff = (np.interp(E,res[0], res[1]))
                return eff
# Gas Excitation Cross Section [m2]
def Gas_Excitation(E,species_name):
    with open(CS_dir) as fp:
        processes = parser.parse(fp)
        results= list()
        for i in processes:
            if (i['target']) == species_name and i['kind'] == 'EXCITATION':
                res = [list(z) for z in zip(*i['data'])]
                exc = (np.interp(E,res[0], res[1]))
                results.append(exc)
        return results
# Total Gas Cross Section [m2]
def Sum_Gas_Excitation(E,species_name):
    cross = Gas_Excitation(E,species_name)
    total = np.zeros(len(E))
    for i in range(0,len(cross)):
        total += cross[i]
    return (total)
# Gas Ionization Cross Section [m2]
def Gas_Ionization(E,species_name):
    with open(CS_dir) as fp:
        processes = parser.parse(fp)
        for i in processes:
            while (i['target']) == species_name and i['kind'] == 'IONIZATION':
                res = [list(z) for z in zip(*i['data'])]
                ion = (np.interp(E,res[0], res[1]))
                return ion
# Gas Inelastic Cross Section [m2]            
def Gas_Inelastic(E, species_name):
    return (Sum_Gas_Excitation(E, species_name) + Gas_Ionization(E, species_name))
# Gas Elastic Cross Section [m2]            
def Gas_Elastic(E, species_name):
    return (Gas_Effective(E,species_name) - Gas_Inelastic(E, species_name))
''' ######################## Particle Cross Section ####################### '''
# Elastic Cross Section [m2]: Nanoparticles
def elas(Vp, lam_f, x):
    result = list()
    for i in x:
        if i > 0:
            a =  (2 * np.pi * (rp ** 2) * ((Vp/(2*i)) ** 2))
            b = (((lam_f) ** 2) + (((Vp/(2*i)) ** 2) * (rp ** 2)))
            c = (rp ** 2) * ((((((Vp/(2 * i)) ** 2)))))
            result.append(a * np.log(b/c) + (np.pi * ((rp) ** 2) * (1 + (Vp/(i)))))
        else:
            s = 0.01
            d =  (2 * np.pi * (rp ** 2) * ((Vp/(2*s)) ** 2))
            e = (((lam_f) ** 2) + (((Vp/(2*s)) ** 2) * (rp ** 2)))
            f = (rp ** 2) * ((((((Vp/(2 * s)) ** 2)))))
            result.append(d * np.log(e/f)+ (np.pi * ((rp) ** 2) * (1 + (Vp/(s))))) 
    return result
# Gaussian Potential Distribution described with a mean and a standard deviation
def vp_dist(mean, x, e_temperature, n_e, n_i):
    result = list()
    v = n_e / n_i
    t = (e_temperature * 11604.52500617) / Ti
    nu = (co.m_e / mass_one)
    charge = (4 * np.pi * co.epsilon_0 * rp * mean) / co.e
    be = np.power(co.e,2) /(4*np.pi*co.epsilon_0*rp*co.k*e_temperature*11604.52500617)
    sigma2q  = ((1/be) *(1 - t*be*charge)) / (t + 0.899*v*np.sqrt(t/nu)*sp.special.erfc(-0.494*be*charge))
    sigma = (np.sqrt(sigma2q) *co.e)/ (4 * np.pi * co.epsilon_0 * rp)
    vp_dist = (1/(np.sqrt(2*np.pi*sigma)))*np.exp(-0.5*(((x + mean)/sigma)**2))
    vp_area = np.trapz(vp_dist,x)
    vp_norm = vp_dist / vp_area
    return (mean, 3*0.2)


#print('sigma', vp_dist(-0.9809265081394956,E,5.519222065628645,1.38e+15,1.2278919316041084e+16)[1]/3)
# Attachment (Inelastic) Cross Section [m2]: Nanoparticles
def attach(Vp, x):
    result = list()
    for i in x:
        if i > -Vp:
            result.append((np.pi * ((rp) ** 2) * (1 + (Vp/(i)))))
        else:
            result.append(0)
    return result  
# Sum of attachment cross section with respect to poteential distribution 
def attach_dist(mean, x, e_temperature, n_e, n_i):
    res = []
    vp = vp_dist(mean, x, e_temperature, n_e, n_i)
    num = int(round(2*vp[1]*100))
    vp_range = np.linspace(vp[0] - vp[1], vp[0] + vp[1],num)
    for k in vp_range:
        res.append(attach(k,x))
    return np.sum(res, axis=0)
# Total Kinetic energy of ejected electrons: removal energy + particle Potential 
def removal_energy(charge):
    res = W_bulk/(rp*1E9*10)
    electro = ((2*np.absolute(charge)-1) * co.e) / (8 * np.pi * co.epsilon_0 * rp)
    rem = W_bulk + res - electro
    if rem > 0:
        return rem
    else:
        return 0
# Detachment Cross Sections [m2]: Nanoparticles
def detach(Vp, x, prop_factor):
    hump_eng = (-E_fermi) + Vp
    result = list()
    for i in x:
            result.append((1 * prop_factor * (1/((2 * np.pi * 0.001)) ** 0.5) *
                           np.exp(-((i + 1 *(hump_eng)) ** 2)/(0.2**2))))
    return result       


#plt.legend(['Attachment', 'Detachment'])

E = np.linspace(0,30,2999)
path_dir = r'/Users/scienceman/Desktop/Project_Titan/Falcon/Dusty_Day_One/blackbox_photoemission/'

T = pd.read_csv(path_dir + 'EEPF_10Pa_dist_photo.txt')
EEPF_1 = np.array(T.iloc[:,1]).astype(float)
  
fig, ax1 = plt.subplots()
fig.set_size_inches(10, 10)

#poteentialsss = -0.9376544468780817
char_conv = ((4*np.pi*co.epsilon_0*rp* 0.9376544468780817)/co.e)


print(char_conv)
#ax1.rcParams.update({'font.size' : 13})

color = 'tab:red'
ax1.set_xlabel('Energy (eV)', fontsize=13)
ax1.set_ylabel('Cross section (m2)', fontsize=13)
ax1.tick_params(which='minor', direction='in')
ax1.tick_params(direction= 'in')
#ax1.plot(p,temp, color='r', linewidth=3, marker='o')
ax1.plot(E, attach_dist(-1.8924240176024114,E,4,1e15,1e16), linewidth=2, color='b')
ax1.plot(E, detach(-1.8924240176024114, E,0.57E-13 ), linewidth=2, color='r')
ax1.legend(['Attachment', 'Detachment'])
ax1.text(1,7e-13, r'SD: 0.2 eV - 75 mTorr', fontsize=13)


#ax1.set_ylim(1000,2000)
#plt.xticks(np.arange(30,120, step=10))
#plt.xlim(28,110)
#ax1.tick_params(axis='y', labelcolor=color)
plt.grid(True)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#ax1.text(30,1110, r'$(c)$',fontsize=13)
color = 'tab:blue'
ax2.set_ylabel('$f_{o}$$(eV^{-3/2})$', color='g', fontsize=13)  # we already handled the x-label with ax1
ax2.tick_params(which='minor', direction='in')
ax2.tick_params(direction= 'in')
ax2.plot(E, EEPF_1, linewidth=3, color= 'g')
ax2.set_yscale('log')
ax2.set_ylim(0.01,0.06)
ax2.set_xlim(0,16)
#ax2.plot(p,thii, color='lime', linewidth=3, marker= 'o')
#ax2.plot(p,photoii, color='gold', linewidth=3, marker= 'o')
#ax2.plot(p,qii, color='b', linewidth=3, marker= 'o')
#ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()













