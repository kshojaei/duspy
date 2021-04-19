from __future__ import division
import numpy as np 
import pandas as pd
from bolos import solver, grid, parser, process, target
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as co
import scipy as sp
import miepython
'''############# Input Values ############## '''
ne_guess = 3.872e15                # Electron Density [1/m3] (Initial Guess)
Vp_guess = -0.19524391310775338# Particle Potential (eV) (Initial Guess)
ni_guess = 1.1197754158125372e+16 # Ion Density [1/m3]
EN_range = np.linspace(0,1,int(1E1))# Range of iterable Reduced Electric Field 
ni_range = np.linspace(0, 1E17,int(1E3))
amp_range = np.linspace(0,-1E-17,int(1E3))
rp = 16E-9                          # Particle Radius [m]
RR = 0.1                            # Maximum Discharge Dimension [m]
npp = 5E15                          # Particle Density [1/m3]
Te_guess = 5.6                      # Electron Temperature [eV]
EN_guess = 760.7           # Reduced Electric Field [Td]
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
E = np.linspace(0,30,3000)          # Electron Energy Range [eV]
bb_guess = 11.65E-15               # Amplitude of detachment cross section
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
sd = 0.05                           # Standard deviation of potential distribution
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
g_escape_1 = 2E-3/5               # Escape Factor (3P1) Ex Argon lmda_1: 106.7E-9 nm
g_escape_2 = 2.5E-3/5              # Escape Factor (1P1) Ex Argon lmda_2: 104.8E-9 nm
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
# Potential Distribution
def weight(mean, s,Vp):
    one = 1/(np.sqrt(2*np.pi)*s)
    two = ((Vp - mean)/s)**2
    return one * np.exp(-0.5 * two)
# Attachment (Inelastic) Cross Section [m2]: Nanoparticles for one potential
def attach(Vp, x):
    result = list()
    for i in x:
        if i > -Vp:
            result.append((np.pi * ((rp) ** 2) * (1 + (Vp/(i)))))
        else:
            result.append(0)
    return result 
# Attachment Cross Section with respect to potential distribution 
def attach_dist(mean,x,s):
    result = list()
    vp_range = np.linspace(-8,-0.01,100)
    for i in vp_range:
        result.append(np.multiply(weight(mean, s, i),attach(i,x)))
    return sum(result) * (7.99/100)
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
def detach(Vp, x, prop_factor,s):
    hump_eng = (-E_fermi) + Vp
    result = list()
    for i in x:
            result.append((-1 * prop_factor * (1/((2 * np.pi * 0.001)) ** 0.5) *
                           np.exp(-((i + 1 *(hump_eng)) ** 2)/(s**2))))
    return result 
# Electron Energy Relaxation Length [m]
def EERL(E, species_one, species_two, f, e_temperature, e_density, i_density, floating_potential):
    pn = (p/ (co.R * T_gas)) * co.Avogadro
    # Electron Debye Length [m]
    lamd_e = np.sqrt((co.epsilon_0 * co.k * e_temperature * 11604.52500617 )/
                (co.e * co.e * e_density))
    # Ion Debye Length [m]
    lamd_i = np.sqrt((co.epsilon_0 * co.k * Ti)/(co.e * co.e * i_density))
    # Particle Temperature [K]
    T_particle = 300 
    # Particle Charge
    kk = (i_density - e_density)/(npp)
    # Particle Debye Length [m]
    lamd_p = np.sqrt((co.epsilon_0 * co.k * T_particle)/(kk * co.e * co.e * npp))
    # Linearized Debye Length [m]
    lamd_f = ((1/lamd_e)+(1/lamd_i)+(1/lamd_p))**(-1)
    # Species One: Inverse Mean Free Path [1/m] (Inelastic)  
    d_one_in = (f * pn) * ((Gas_Elastic(E, species_one) *
                          (2*co.m_e/mass_one)) +  Gas_Inelastic(E, species_one))
    # Species Two: Inverse Mean Free Path [1/m] (Inelastic)
    d_two_in = ((1-f) * pn) * ((Gas_Elastic(E, species_two) * 
                (2*co.m_e/mass_two)) +  Gas_Inelastic(E, species_two))
    # Particle: Inverse Mean Free Path [1/m] (Inelastic)
    d_p_in = npp * (np.multiply((2 * co.m_e / m_particle), elas(floating_potential,
                                lamd_f,E)) + attach_dist(floating_potential, E,sd))
    # Species One: Inverse Mean Free Path [1/m] (Momentum Transfer)  
    d_one_el = (f * pn) * Gas_Elastic(E, species_one)
    # Species Two: Inverse Mean Free Path [1/m] (Momentum Transfer)  
    d_two_el = ((1- f) * pn) * Gas_Elastic(E, species_two) 
    # Particle: Inverse Mean Free Path [1/m] (Momentum Transfer)
    d_p_el = np.multiply(npp, elas(floating_potential,lamd_f,E))
    # Total Mean Free Path for Momentum Transfer [m] With Particle
    mfp_mp = 1 / (d_one_el + d_two_el + d_p_el)
    # Total Mean Free Path for all Collisional Energy Loss Processes [m] With Particles
    mfp_inelasticp = 1 / (d_one_in + d_two_in + d_p_in)
    # Total Mean Free Path for Momentum Transfer [m] Without Particle
    mfp_m = 1 / (d_one_el + d_two_el)
    # Total Mean Free Path for all Collisional Energy Loss Processes [m] Without Particles
    mfp_inelastic = 1 / (d_one_in + d_two_in)
    # Electron Energy Relaxation Length With Particles [m]
    rlp = np.sqrt((mfp_mp * mfp_inelasticp)/3)    
    # Electron Energy Relaxation Length Without Particles [m]
    rl = np.sqrt((mfp_m * mfp_inelastic)/3) 
    return rlp, rl
''' ########################## Laser ###############################'''
class Laser:
    def __init__(self, laser_power, laser_wavelength, laser_spot_radius, particle_radius, mat_dir):
        self.laser_power = laser_power
        self.laser_wavelength = laser_wavelength
        self.laser_spot_radius = laser_spot_radius
        self.particle_radius = particle_radius
        self.mat_dir = mat_dir
    # Spot Size Calculator
    def spot_size(self, focal_length):
        return (1.27 * focal_length * self.laser_wavelength * (mtwo / (2*laser_spot_radius)))
    # Absorption, Scattering, and Exrinction Cross Section of a Particle [m2]
    def abs_scatt_ext(self):
        # import from material properties from DataBase
        X = pd.read_csv(self.mat_dir, delimiter = ',')
        # Wavelength (nm)
        ref_lam = np.array(X.iloc[:,0]).astype(float)
        # real part of refractive index 
        ref_n =  np.array(X.iloc[:,1]).astype(float)
        # imaginary part of refractive index
        ref_k =  np.array(X.iloc[:,2]).astype(float)
        # Interpolate the real part of refractive index
        n_interp = np.interp((self.laser_wavelength*1E6),ref_lam, ref_n)
        # Interpolate the imaginary part of refractive index
        k_interp = np.interp((self.laser_wavelength*1E6),ref_lam, ref_k)
        # refractive index of material
        ref_m = n_interp - 1.0j*k_interp
        # Particle Radius [microns]
        ref_rp = self.particle_radius * (1E9/1E3)
        # Size Parameter 
        ref_x = (2 * np.pi * ref_rp) / (self.laser_wavelength*1E6)
        qext, qsca, qback, g = miepython.mie(ref_m,ref_x)
        # Absorption Efficiency
        qabs = qext - qsca
        # Absorption Cross Section (m2)
        qqabs = qabs * np.pi * (self.particle_radius**2)
        # Scattering Cross Section 
        qqsca = qsca * np.pi * (self.particle_radius**2)
        # Extinction Cross Section 
        qqext = qqsca + qqabs
        return qqabs, qqsca, qqext
    # Power absorbed by a particle from a CW laser (W)
    def cw_laser_absoprtion(self):
        # CW Laser Irradiance of a Gaussian profile
        pow_density = self.laser_power / (np.pi * (self.laser_spot_radius**2)/2)
        return (pow_density * self.abs_scatt_ext()[0])
    
    def pulsed_laser_absorption(self, laser_energy, repetition_rate):
        max_pow_density = ((laser_energy*repetition_rate) / 
                           (np.pi * (self.laser_spot_radius**2)/2))
        return (max_pow_density * self.abs_scatt_ext()[0])
''' ############### Argon Resonance Radiation and Photoemission ############ '''
class Resonance:
    def __init__(self, electron_density, ion_density, EEPF, charge, species_name):
        self.electron_density = electron_density
        self.ion_density = ion_density
        self.EEPF = EEPF
        self.species_name = species_name
        self.charge = charge
    # Normalized Electron Energy Distribution Function
    def EEDF_norm(self):
        return ((self.EEPF * np.sqrt(E)) / np.trapz((self.EEPF * np.sqrt(E)),E))  
    def Quantum_Probs(self):
        vo_Ar = np.sqrt(2 * co.R * T_gas / M_Ar)  # Average Argon Velocity [m/s]
        pn = (p/ (co.R * T_gas)) * co.Avogadro    # Gas Density  [1/m3]
        ##################################
        # Resonance Line: The line of longest wavelength associated with a transition
        # between the ground state and an excited state
        lmda_1 = 106.7E-9                # Resonance line (3P1) excited argon [m]
        lmda_2 = 104.8E-9                # Resonance line (1P1) excited argon [m]
        gama_1 = 1.32E8                  # Transition Probability (3P1) [1/s]
        gama_2 = 5.32E8                  # Transition Probability (1P1) [1/s]
        vo_1 = co.c / lmda_1             # frequency at which k(v) is max (3P1) Ex Argon
        vo_2 = co.c / lmda_2             # frequency at which k(v) is max (1P1) Ex Argon
        tau_1 = 1 / gama_1               # Lifetime of the Ex Argon (3P1) [s]
        tau_2 = 1 / gama_2               # Lifetime of the Ex Argon (1P1) [s]
        g_1 = 3                          # Degeneracy (3P1) Ex Argon (2L+1)(2S+1)
        g_2 = 3                          # Degeneracy (1P1) Ex Argon 
        g_0 = 1                          # Degeneracy Argon Normal State (1S0)
        # Spectroscopic Notation: (2S+1)LJ
        k_01 = (((lmda_1 **2)* pn)/(8*np.pi))*(g_1/g_0)*(1/tau_1)
        k_02 = (((lmda_2 **2)* pn)/(8*np.pi))*(g_2/g_0)*(1/tau_2)
        f_1 = 0.0675                     # Oscillator Strength for 106.7 nm (3P1) Ex Argon
        f_2 = 0.2629                     # Oscillator Strength for 104.8-9 nm (1P1) Ex Argon
        # Average collision rate for the Ex Argon (3P1) [1/s]
        gama_c_1 = (4 * co.e * co.e * f_1 * lmda_1 * pn)/(3 * mass_one * co.c)
        # Average collision rate for the Ex Argon (1P1) [1/s]
        gama_c_2 = (4 * co.e * co.e * f_2 * lmda_2 * pn)/(3 * mass_one * co.c)
        a_1  = ((gama_1 + gama_c_1) * lmda_1) / (4 * np.pi * vo_Ar)
        a_2  = ((gama_2 + gama_c_2) * lmda_2) / (4 * np.pi * vo_Ar)
        # Impressed Frequency [1/s]
        v_1 = np.linspace(vo_1 -1E11, vo_1 + 1E11,3000)
        v_2 = np.linspace(vo_2 -1E11, vo_2 + 1E11,3000)
        x_1 = (v_1 - vo_1)*(lmda_1 / vo_Ar)
        x_2 = (v_2 - vo_2)*(lmda_2 / vo_Ar)
        def Voig_1(u,y):
            top =  np.exp(-(y**2))
            bottom = a_1**2 + (u - y)**2
            return top/bottom 
        def Voig_2(u,y):
            top =  np.exp(-(y**2))
            bottom = a_2**2 + (u - y)**2
            return top/bottom 
        # Voigt Profile [3P1] Ex Argon
        def Voigt_1(u):
            res = np.zeros_like(u)
            for i,val in enumerate(u):
                y,err = sp.integrate.quad(Voig_1, -np.inf, np.inf, args=(val))
                res[i] = y
            return (a_1/np.pi) * res
        # Voigt Profile [1p1] Ex Argon
        def Voigt_2(u):
            res = np.zeros_like(u)
            for i,val in enumerate(u):
                y,err = sp.integrate.quad(Voig_2, -np.inf, np.inf, args=(val))
                res[i] = y 
            return (a_2/np.pi) * res
        Voigt_1_area = np.trapz(Voigt_1(x_1),v_1)
        Voigt_2_area = np.trapz(Voigt_2(x_2),v_2)
        # Normalized Voigt_1 Function [1P1] Ex Argon
        Voigt_1_norm = Voigt_1(x_1) / Voigt_1_area
        # Normalized Voigt_2 Function [3P1] Ex Argon
        Voigt_2_norm = Voigt_2(x_2) / Voigt_2_area
        # Absorption Coefficients accounting for natural, Doppler, and 
        # pressure broadenings (3P1) Ex Argon (KP(v))
        kv_1 = k_01 * Voigt_1_norm   
        # Absorption Coefficients accounting for natural, Doppler, and 
        # pressure broadenings (1P1) Ex Argon (KP(v))
        kv_2 = k_02 * Voigt_2_norm  
        # Probability of a quantum penentrating a distance rho in the gas 
        # before being absorbed for constant rho 
        def T_11(x):
            result = list()
            for i in x: 
                if i >= 0:
                    result.append(np.trapz((Voigt_1_norm * np.exp(-kv_1 * i)),v_1))
                else:
                    result.append(0)
            return result
        # Probability of quantum penetrating a distance rho in a gas
        # before being absorbed for series of rho
        def T_22(x):
            result = list()
            for i in x: 
                if i >= 0:
                    result.append(np.trapz((Voigt_2_norm * np.exp(-kv_2 * i)),v_2))
                else:
                    result.append(0)
            return result
        radial = np.linspace(0, RR, 3000)
        return (np.trapz(T_11(radial),radial)) , (np.trapz(T_22(radial),radial))
    # photoemission frequency by UV photodetachment of Ar (3P1) & (1P1) (1/s)
    def photoemission_freq(self):
        with open(CS_dir) as fp:
            processes = parser.parse(fp)
            for i in processes:
                while (i['target']) == self.species_name and i['kind'] == 'EXCITATION':
                    res = [list(z) for z in zip(*i['data'])]
                    cs = (np.interp(E,res[0], res[1]))
                    pn = (p/ (co.R * T_gas)) * co.Avogadro
                    vex = pn * cs * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
                    vexc = np.trapz(vex,E)
                    #Population Density of Ex Argon (3P1) & (1P1)
                    n_1 = (self.electron_density * vexc) / (g_escape_1 * gama_1) 
                    n_2 = (self.electron_density * vexc) / (g_escape_2 * gama_2)
                    # Charging Frequency by resonance UV photodetachment
                    freq_g_1 = Q * np.pi * (rp**2) * n_1 * gama_1 * self.Quantum_Probs()[0]
                    #print('Quantum Prob 1', self.Quantum_Probs()[0])
                    
                    #freq_g_1 = Q * np.pi * (rp**2) * n_1 * gama_1 * 0.0001135311880446913
                    freq_g_2 = Q * np.pi * (rp**2) * n_2 * gama_2 * self.Quantum_Probs()[1]
                    #print('Quantum Prob 2', self.Quantum_Probs()[1])
                    #freq_g_2 = Q * np.pi * (rp**2) * n_2 * gama_2 * 3.78660267195606e-05
                    return (freq_g_1 + freq_g_2)
        # Charging Frequency Due to Quenching Process          
    def quench_freq(self):
        with open(CS_dir) as fp:
            processes = parser.parse(fp)
            for i in processes:
                while (i['target']) == self.species_name and i['kind'] == 'EXCITATION':
                    res = [list(z) for z in zip(*i['data'])]
                    cs = (np.interp(E,res[0], res[1]))
                    pn = (p/ (co.R * T_gas)) * co.Avogadro
                    vex = pn * cs * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
                    vexc = np.trapz(vex,E)
                    s_1 = (self.electron_density * vexc) / (g_escape_1 * gama_1) 
                    s_2 = (self.electron_density * vexc) / (g_escape_2 * gama_2)
                    return ((s_1 + s_2) * np.sqrt((co.k * T_gas)/mass_one) * np.pi * (rp ** 2)) 
    # Ion Frequency [1/s]
    def freq_i(self):
        S = 4 * np.pi * (rp**2)
        Vp = (-(self.charge * co.e) / ( 4 * np.pi * co.epsilon_0 * rp))
        A = S * (self.ion_density) * np.sqrt((co.k * T_gas)/(2 * np.pi * mass_dom))
        nu = A * (1 - ((co.e * Vp)/(co.k * T_gas)))
        return nu
    def photoion(self):
        return (self.photoemission_freq() / self.freq_i())
    def quenchion(self):
        return (self.quench_freq() / self.freq_i())
''' ########################## Heat Balance ############################### '''
class Heat:
    def __init__(self, charge, electron_temperature, EEPF, ion_density, X_edii, X_edea):
        self.charge = charge   
        self.X_edii = X_edii
        self.X_edea = X_edea
        self.electron_temperature = electron_temperature
        self.EEPF = EEPF
        self.ion_density = ion_density
        ''' HEAT LOSS MECHANISMS '''
    # Normalized Electron Energy Distribution Function
    def EEDF_norm(self):
        return ((self.EEPF * np.sqrt(E)) / np.trapz((self.EEPF * np.sqrt(E)),E))  
    # Condution [J/s]
    def conduction(self,Tp):
        surf = np.pi * (rp ** 2)
        KE = 1.5 * co.k * (Tp - T_gas)
        pn = (p/ (co.R * T_gas)) * co.Avogadro
        one = surf * KE * f_maxmin * pn * np.sqrt((8 * co.k * T_gas)/(np.pi * mass_one))
        two = surf * KE * (1-f_maxmin) * pn * np.sqrt((8 * co.k * T_gas)/(np.pi * mass_two))
        return (one + two)
    # Radiation Loss [J/s]: Gold
    def rad_gold(self,Tp):
        # import from DataBase: Gold
        T = pd.read_csv(gold_dir, delimiter = ',')
        # Wavelength (nm)
        lamg = np.array(T.iloc[:,0]).astype(float)*1E-6 
        # ng 
        ng =  np.array(T.iloc[:,1]).astype(float)
        # kg
        kg =  np.array(T.iloc[:,2]).astype(float)
        Eg = (6*ng*kg)/((((ng**2) - (kg**2) + 2 )**2) + 4*((ng*kg)**2))
        emis_g =(4*np.pi*(2*rp)*Eg) / (lamg)
        result = list()
        for i in Tp:
            if i > 0:
                rgold = emis_g * np.pi*((2*rp)**2)*(2 * np.pi* co.h* co.c* co.c)/((lamg**5)*
                                        (np.exp((co.h*co.c)/(lamg*co.k*i)) - 1))            
                result.append((np.trapz(rgold,lamg)))
            else:
                result.append(0)
        return result
    # Radiation Loss [J/s]: Copper
    def rad_copper(self,Tp):
        # import from DataBase: Copper
        X = pd.read_csv(copper_dir, delimiter = ',')
        # Wavelength (nm)
        lamgc = np.array(X.iloc[:,0]).astype(float)*1E-6 
        # ng 
        ngc =  np.array(X.iloc[:,1]).astype(float)
        # kg
        kgc =  np.array(X.iloc[:,2]).astype(float)
        Egc = (6*ngc*kgc)/((((ngc**2) - (kgc**2) + 2 )**2) + 4*((ngc*kgc)**2))
        emis_gc =(4*np.pi*(2*rp)*Egc) / (lamgc)
        #result = list()
        #for i in Tp:
        #    if i > 0:
        rcopper = emis_gc * np.pi*((2*rp)**2)*(2 * np.pi* co.h* co.c* co.c)/((lamgc**5)*
                                                                            (np.exp((co.h*co.c)/(lamgc*co.k*Tp)) - 1))            
#                result.append((np.trapz(rcopper,lamgc)))
#            else:
#                result.append(0)
        return (np.trapz(rcopper,lamgc))
    # Radiation Loss [J/s] Carbon
    def rad_carbon(self, Tp):
        lamda = np.linspace(10E-9,100000E-9,3000)   
        nm = 2.213 + 9.551E3 * lamda
        km = 0.7528 + 1.265E4 * lamda
        Em = (6 * nm * km)/((((nm ** 2) - (km ** 2) + 2 )**2) + 4*((nm * km) ** 2))
        emis = (4 * np.pi * (2 * rp) * Em) / (lamda)
        result = list()
        for i in Tp:
            if i > 0:
                radd = emis*np.pi*((2*rp)**2)*(2*np.pi*co.h*co.c*co.c)/((lamda**5)*
                                   (np.exp((co.h*co.c)/(lamda*co.k*i)) - 1))            
                result.append((np.trapz(radd,lamda)))
            else:
                result.append(0)
        return result
    # Removal [J/s]
    def rem_func(self):
        res = W_bulk/(rp*1E9*10)
        electro = ((2*np.absolute(self.charge)-1) * co.e) / (8 * np.pi * co.epsilon_0 * rp)
        rem = W_bulk + res - electro 
        if rem > 0:
            return rem
        else:
            return 0
        '''
    def freq_i(self, mean,x, e_temperature, n_e, n_i):
        res = []
        sa = 4 * np.pi * (rp ** 2)
        a_coeff = sa * n_i * (((co.k*Ti)/(2* np.pi * mass_one))** 0.5)
        vp = vp_dist(mean, x, e_temperature, n_e, n_i)
        num = int(round(2*vp[1]*100))
        vp_range = np.linspace(vp[0] - vp[1], vp[0] + vp[1],num)
        for k in vp_range:
            res.append(sa * a_coeff * np.exp(-((co.e * k)/(co.k * Ti))))
        return np.sum(res)
    '''               
    # Ion Frequency [1/s]
    def freq_i(self):
        S = 4 * np.pi * (rp**2)
        Vp = (-(self.charge * co.e) / ( 4 * np.pi * co.epsilon_0 * rp))
        A = S * (self.ion_density) * np.sqrt((co.k * T_gas)/(2 * np.pi * mass_dom))
        nu = A * (1 - ((co.e * Vp)/(co.k * T_gas)))
        return nu
    # Thermionic Emission [J/s]
    def Thermionic(self):
        rem_nano = co.e * self.rem_func()
        #prop = 4 * co.m_e * ((2 * rp * np.pi * co.k * Tp)**2)
        #return (prop/((co.h)**3)) * rem_nano * np.exp(- rem_nano / (co.k * Tp))   
        return (rem_nano * self.X_edii * self.freq_i())
    # Thermionic Frequency
    def thermion(self, Tp):
        rem_nano = co.e * self.rem_func()
        Coeff = 4 * co.m_e * ((2*rp*np.pi*co.k*Tp)**2)
        thermf = ((Coeff/((co.h)**3)) * 1 * np.exp(- rem_nano / (co.k * Tp)))
        return (thermf / self.freq_i())
    def detachion(self):
        return (self.X_edii)
    # Dissociative  and Associative Sequence (H3+ -> H2 + H+ & H+ + H+ -> H2) Heat Gain [J]
    def AD(self):
        # Ionization Potential of Hydrogen (eV)
        E_inz = 13.6
        # Bond Energy of the Hydrogen molecule
        E_diss = 4.52
        return (self.freq_i() * co.e * (E_inz + (E_diss/2)))
    # Electron Kinetic Power [J/s]
    def EKE(self):
        potential = (co.e * abs(self.charge))/(4 * np.pi * co.epsilon_0 * rp)
        h_index = (np.argwhere((E)> potential))[0][0]
        f_area = np.trapz(self.EEDF_norm()[h_index:], E[h_index:])
        normf = (self.EEDF_norm()[h_index:])/f_area
        KE_e = np.trapz((E[h_index:]*normf), E[h_index:])
        #return ((self.X_edii / self.X_edea) * self.freq_i() * KE_e) * co.e
        #velocity = np.sqrt((2*E[h_index])/co.m_e)
        #up = 0.5 * co.m_e * (velocity**3) * normf
        #bo = velocity * normf
        #up_a = np.trapz(up, velocity)
        #bo_a = np.trapz(bo, velocity)
        #return ((self.X_edii / self.X_edea) * self.freq_i() * (up_a/bo_a) * co.e)
        return ((self.X_edii / self.X_edea) * self.freq_i() * (2 * co.k *
                self.electron_temperature * 11604.52500617))
    # Ion Kinetic Power [J/s]
    def IKE(self):
        Vp = np.abs((- self.charge * co.e) / ( 4 * np.pi * co.epsilon_0 * rp))
        return (self.freq_i()*((0.5*co.k * self.electron_temperature * 11604.52500617)+ (Vp*co.e)))
    # Atomic Hyrogen-Particle Collision Heating [J/s]
    def AtomicH(self):
        Sa = 4 * np.pi * (rp**2)
        # bond energy of the hydrogen molecule (eV)
        hbond_energy = 4.51
        ion_den = self.ion_density * 1E2
        # partilce_radical hydrogen collision frequency (1/s)
        hfreq = 0.25 * ion_den * Sa * np.sqrt((8*co.k*T_gas)/(np.pi* mass_H))
        return (hfreq * hbond_energy * co.e) 
    # Metastable Argon 3P2 Impaction Energy Rate [J/s]
    def E_3p2(self):
        n_3p2 = (2E10/1E11) * self.ion_density
        surf = (rp**2) * np.pi 
        v_avg = np.sqrt((8 * co.k * T_gas)/(np.pi * mass_one))
        Em3p2 = 11.548
        return (n_3p2 * surf * v_avg * co.e * Em3p2)
    # Metastable Argon 3P0 Impaction Energy Rate [J/s]
    def E_3p0(self):
        n_3p0 = (3E9/3E11) * self.ion_density
        surf = (rp**2) * np.pi 
        v_avg = np.sqrt((8 * co.k * T_gas)/(np.pi * mass_one))
        Em3p0 = 11.723
        return (n_3p0 * surf * v_avg * co.e * Em3p0)
    # Laser Absorption [J/s]
    def Laser_Heating(self):
        lheat = Laser(laser_pow, laser_wavel, laser_spot_radius, rp, copper_dir).cw_laser_absoprtion()
        return lheat
    # Total Heat Loss: Copper [J/s]
    def Total_Loss_Copper(self, Tp):
        one = self.rad_copper(Tp)
        two = self.Thermionic()
        three = self.conduction(Tp)
        return (one + two + three)
    # Total Heat Loss: Carbon [J/s]
    # Total Heat Loss: Gold [J/s]
    def Total_Loss_Gold(self, Tp):
        one = self.rad_gold(Tp)
        two = self.Thermionic()
        three = self.conduction(Tp)
        return (one + two + three)
    # Total Heat Loss: Carbon [J/s]
    def Total_Loss_Carbon(self, Tp):
        one = self.rad_carbon(Tp)
        two = self.Thermionic()
        three = self.conduction(Tp)
        return (one + two + three) 
    # Total Heat Gain [J/s]
    def Total_Gain(self):
        one = self.EKE()
        two = self.IKE()
        three = self.E_3p2()
        four = self.E_3p0()
        five = self.AD()
        six = self.AtomicH()
        seven = self.Laser_Heating()
        return (one + two + three + four + five + six + seven)
    # Net Heat Rate Carbon [J/s]
    def Energy_Balance_Carbon(self, Tp):
        one = self.Total_Loss_Carbon(Tp)
        two = self.Total_Gain()
        return (two - one)
    # Net Heat Rate Gold [J/s]
    def Energy_Balance_Gold(self, Tp):
        one = self.Total_Loss_Gold(Tp)
        two = self.Total_Gain()
        return (two - one)    
        # Net Heat Rate Gold [J/s]
    def Energy_Balance_Copper(self, Tp):
        one = self.Total_Loss_Copper(Tp)
        two = self.Total_Gain()
        return (two - one)
    #Solving for Particle temperature [K]: Carbon
    def Temp_Solver_Carbon(self):
        return fsolve(self.Energy_Balance_Carbon, 500) 
        #Solving for Particle temperature [K]: Carbon
    def Temp_Solver_Gold(self):
        return fsolve(self.Energy_Balance_Gold, 500)   
    def Temp_Solver_Copper(self):
        return fsolve(self.Energy_Balance_Copper, 500) 
''' ########################## Power Balance ############################### '''
class Power_Balance:
    def __init__(self, EEPF, electron_density, ion_density, Particle_Potential, 
                 electron_temperature, species_name):
        self.electron_density = electron_density 
        self.Particle_Potential = Particle_Potential
        self.EEPF = EEPF
        self.electron_temperature = electron_temperature
        self.ion_density = ion_density 
        self.species_name = species_name
    # Normalized Electron Energy Distribution Function
    def EEDF_norm(self):
        return ((self.EEPF * np.sqrt(E)) / np.trapz((self.EEPF * np.sqrt(E)),E))  
    def EEPF_norm(self):
        return ((self.EEPF) / np.trapz((self.EEPF),E))         
    # Electron-Neutral Momentum Transfer
    def L_momentum(self):
        V = (R ** 2) * np.pi * height
        Kin = 1.5 * co.k * ((self.electron_temperature * 11604.52500617)- T_gas)  
        pn = (p/ (co.R * T_gas)) * co.Avogadro
        vm = pn * Gas_Elastic(E,self.species_name) * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
        vmm = np.trapz(vm, E)
        return (2 * self.electron_density * V * Kin * vmm * (co.m_e / mass_one))
    # Energy loss to ionization events in plasma discharge
    def L_ionization(self):
        V = (R ** 2) * np.pi * height    
        with open(CS_dir) as fp:
            processes = parser.parse(fp)
            for i in processes:
                while (i['target']) == self.species_name and i['kind'] == 'IONIZATION':
                    res = [list(z) for z in zip(*i['data'])]
                    cs = (np.interp(E,res[0], res[1]))
                    pn = (p/ (co.R * T_gas)) * co.Avogadro
                    vion = pn * cs * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
                    vionc = np.trapz(vion, E)
                    return (2 * self.electron_density * vionc * i['threshold'] * V * co.e) 
    # Energy loss to excitation events in plasma discharge        
    def L_excitation(self):
        V = (R ** 2) * np.pi * height 
        results = list()
        with open(CS_dir) as fp:
            processes = parser.parse(fp)
            for i in processes:
                if (i['target']) == self.species_name and i['kind'] == 'EXCITATION':
                    res = [list(z) for z in zip(*i['data'])]
                    cs = (np.interp(E,res[0], res[1]))
                    pn = (p/ (co.R * T_gas)) * co.Avogadro
                    vex = pn * cs * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
                    vexc = np.trapz(vex,E)
                    lexx = self.electron_density * vexc * i['threshold'] * V * co.e 
                    results.append(lexx)
            return np.sum(results)
    # Energy loss through diffusion to the electrodes
    def L_ion_acc(self):
        #pn = (p/ (co.R * T_gas)) * co.Avogadro
        #vm = pn * Gas_Elastic(E,self.species_name) * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()   
        #vmm = np.trapz(vm, E)
        #mo_e = co.e / (co.m_e * vmm)
        #D_e = (co.k * self.electron_temperature * 11604.52500617) / (co.m_e * vmm)
        lambda_i = 1E-2/(330 * (0.0075006168 * p))
        D_i = ((3 * np.pi)/(16 * np.sqrt(2))) * np.sqrt((8 * co.k * Ti)/(np.pi * mass_one)) * lambda_i
        #mo_i = (co.e * D_i)/(co.k * Ti)
        #D_a = (mo_i * D_e + mo_e * D_i)/(mo_e + mo_i)
        D_a = D_i * ((11604.52500617* self.electron_temperature)/ Ti)
        AA = 0*(2 * np.pi * R * height) + (2* np.pi * (R**2))
        return self.ion_density * (D_a/(height/2)) * AA * (Vsh/2) * co.e
    # Energy loss to the nanoparticles through the collection of electrons
    def L_ecollection(self):
        V = (R ** 2) * np.pi * height
        vel = ((1 + (self.Particle_Potential/(E+ 1E-40))) * np.sqrt((2 * co.e)/
                (co.m_e)) * self.EEPF_norm()* E)
        Koeff = (npp * V * co.e) * (np.pi * (rp**2) * self.electron_density)*E
        PR = (np.argwhere((Koeff*vel)>0))[0][0]
        return np.trapz((Koeff * vel)[PR:], E[PR:])
    def Power_Loss(self):
        one = self.L_momentum()
        two = self.L_ionization()
        three = self.L_excitation()
        four = self.L_ion_acc()
        five =  self.L_ecollection()
        return one + two + three + four + five
    def Net_Power(self):
        return (self.Power_Loss() - P_in)/ P_in
''' ######################### Particle Balance ############################ '''
class Particle_Balance:
    def __init__(self, EEPF, species_name):
        self.EEPF = EEPF
        self.species_name = species_name
    # Normaized Electron Energy Distribution Function (eV-1 vs. eV)
    def EEDF_norm(self):
        return ((self.EEPF * np.sqrt(E)) / np.trapz((self.EEPF * np.sqrt(E)),E)) 
    # Rate of Ion Generation  (Particle/s)
    def Ion_Generation(self, electron_density):
        V = (R ** 2) * np.pi * height    
        with open(CS_dir) as fp:
            processes = parser.parse(fp)
            for i in processes:
                while (i['target']) == self.species_name and i['kind'] == 'IONIZATION':
                    res = [list(z) for z in zip(*i['data'])]
                    cs = (np.interp(E,res[0], res[1]))
                    pn = (p/ (co.R * T_gas)) * co.Avogadro
                    vion = pn * cs * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()
                    vionc = np.trapz(vion, E)
                    return (electron_density * vionc * V)
    # Rate of Ion Recombination due to collection via nanoparticls (Particles/s)
    def Ion_Dust_Recomb(self, ion_density, particle_potential):
        V = (R ** 2) * np.pi * height 
        S = 4 * np.pi * (rp**2)
        A = S * (ion_density) * (((co.k * T_gas)/(2 * np.pi * mass_one))** 0.5)
        nu = A * (1 - ((co.e * particle_potential)/(co.k * T_gas)))
        return (npp * nu * V)
    # Rate of Ion Recombination due to diffusion to the chamber wall (Particles/s)
    def Ion_acc_Recomb(self, ion_density, electron_temperature):
        #pn = (p/ (co.R * T_gas)) * co.Avogadro
        #vm = pn * Gas_Elastic(E,self.species_name) * np.sqrt((2*E*co.e)/ co.m_e) * self.EEDF_norm()   
        #vmm = np.trapz(vm, E)
        #mo_e = co.e / (co.m_e * vmm)
        #D_e = (co.k * electron_temperature * 11604.52500617) / (co.m_e * vmm)
        lambda_i = 1E-2/ (330 * (0.0075006168 * p))
        D_i = ((3 * np.pi)/(16 * np.sqrt(2))) * np.sqrt((8 * co.k * Ti)/(np.pi * mass_one)) * lambda_i
        #mo_i = (co.e * D_i)/(co.k * Ti)
        #D_a = (mo_i * D_e + mo_e * D_i)/(mo_e + mo_i)
        D_a = D_i * ((11604.52500617* electron_temperature)/ Ti)
        AA = 0*(2 * np.pi * R * height) + ( 2* np.pi * (R**2))
        return (ion_density * (D_a/(height/2)) * AA)
    # Total Rate of Ion Recombination (Particles/s)
    def Ion_Total_Recomb(self, ion_density, particle_potential, electron_temperature):
        return (self.Ion_Dust_Recomb(ion_density, particle_potential) + self.Ion_acc_Recomb(ion_density,
                electron_temperature))
    # Net Rate of Ion Production (Particles/s)
    def Net_Production_Rate(self, electron_density, ion_density, particle_potential, electron_temperature):
        top = (self.Ion_Generation(electron_density) - self.Ion_Total_Recomb(ion_density,
               particle_potential, electron_temperature))
        bottom = self.Ion_Total_Recomb(ion_density, particle_potential, electron_temperature)
        return (top/bottom)
''' ######################### Boltzmann Solver ############################ '''
class Boltzmann_Solver:
    def __init__(self, EN_power, electron_temperature, ion_density, bb_g, species_name):
        self.species_name = species_name
        self.EN_power = EN_power
        self.ion_density = ion_density
        self.electron_temperature = electron_temperature
        self.bb_g = bb_g
    def EEPF_Vp_ne(self):
        titan = solver.BoltzmannSolver(gs)
        with open('Ar_H2.txt') as fp:
            processes = parser.parse(fp)
        titan.load_collisions(processes)
        for Vpp in np.linspace(4, -15):
            for nee in np.linspace(1E12, 1E18):
                Vp = Vp_guess
                ne = ne_guess
                # Electron Debye Length [m]
                lam_e = np.sqrt((co.epsilon_0 * co.k * self.electron_temperature * 11604.52500617 )/
                                (co.e * co.e * ne))
                # Ion Debye Length [m]
                lam_i = np.sqrt((co.epsilon_0 * co.k * Ti)/(co.e * co.e * self.ion_density))
                # Linearized Debye Length [m]
                particle_temp = 300      #[K]
                pn = (p/ (co.R * T_gas)) * co.Avogadro
                x = npp /(pn + npp)                       # mole fraction
                kk = (self.ion_density - ne)/(npp)
                lam_p = np.sqrt((co.epsilon_0 * co.k * particle_temp)/(kk * co.e * co.e * npp))
                lam_f = ((1/lam_e)+(1/lam_i)+(1/lam_p))**(-1)
                # Loading Nanoparticle Cross Sections
                titan.add_process(kind = "EFFECTIVE", target="Particle", mass_ratio = (co.m_e/(m_particle)), 
                                  data=np.c_[E, elas(Vp,lam_f,E)])
                titan.add_process(kind = "ATTACHMENT", target="Particle", 
                                  data= np.c_[E, attach_dist(Vp,E,sd)])
                titan.add_process(kind= "ATTACHMENT", target="Particle", 
                                  data= np.c_[E, detach(Vp,E,self.bb_g,sd)])
                ###############################################################################
                titan.target['Ar'].density = (f_maxmin) * (1- x)
                titan.target['H2'].density = (1 - f_maxmin) * (1-x)
                titan.target['Particle'].density = x
                titan.kT = 300 * co.k / co.eV
                titan.EN = self.EN_power * solver.TOWNSEND
                titan.init()
                z0 = titan.maxwell(self.electron_temperature)
                z1 = titan.converge(z0, maxn=100, rtol= tolerance)
                Te = np.multiply(0.66666666666666666666666, np.trapz(((E ** 1.5)*z1), E))
                freq_detach = np.sqrt((2 * co.e)/co.m_e) * x * E * z1 * detach(Vp, E, self.bb_g,sd)
                freq_attach = np.sqrt((2 * co.e)/co.m_e) * x * E * z1 * attach_dist(Vp, E,sd)
                #vth  = np.sqrt(co.e * E) * np.sqrt(2/co.m_e)
                area_detach = np.trapz(freq_detach,E)
                area_attach = np.trapz(freq_attach,E)
                #prop_factor = area_detach/(-area_attach)  
                S = 4 * np.pi * (rp ** 2)
                A1 = S * self.ion_density * (((co.k*Ti)/(2* np.pi * mass_one))** 0.5)
                # Ion attachment frequency for the potential distribution
                def ion_freq_dist(mean, n_i,s):
                    res = []
                    sa = 4 * np.pi * (rp ** 2)
                    a_coeff = sa * n_i * (((co.k*Ti)/(2* np.pi * mass_one))** 0.5)
                    vp_range = np.linspace(-8,-0.01,30000)
                    for k in vp_range:
                        ifr = a_coeff * (1-((co.e * k)/(co.k * Ti)))
                        res.append(np.multiply(weight(mean,s,k),ifr))
                    return sum(res) * (7.99/30000)
                def equations(p):
                    m, l = p
                    eqn1 = m - self.ion_density - ((4 * l * np.pi * co.epsilon_0 * rp)/(co.e)) * npp
                    area_e = ((area_attach+area_detach)/x) * m
                    eqn2 = ion_freq_dist(l,self.ion_density, sd) - area_e
                    return (eqn1, eqn2)
                def equations_inv(p):
                    m, l = p
                    eqn1 = m - self.ion_density - ((4 * l * np.pi * co.epsilon_0 * rp)/(co.e)) * npp
                    area_e = ((area_attach+area_detach)/x) * m 
                    eqn2 = ion_freq_dist(l,self.ion_density, sd) - area_e
                    return (eqn1, eqn2)
                if Vp >= 0:
                    m, l = fsolve(equations_inv, (1E15,-2.7))
                else:
                    m, l = fsolve(equations, (1E15,-2.7))
                Vpp = l
                nee = m
                print('')
                print('__________Boltzmann_Solver__________')
                print('')
                print("ne(1/m3): {:.2e}".format(nee))
                print('Vp(eV):',Vpp)
                print('Te(eV):', Te)     
                ''' Save the Normalized EEPF'''
                eqn11 = nee - self.ion_density - ((4 * Vpp * np.pi * co.epsilon_0 * rp)/(co.e)) * npp
                eqn22 = A1 * (1 - ((co.e * Vpp)/(co.k * Ti))) - (((area_attach+area_detach)/x) * nee)
                ii = ion_freq_dist(Vpp,self.ion_density, sd)
                ed = ((-area_detach)/x) * ne
                ea = ((area_attach)/x) * ne
                print(eqn11, eqn22)
                print('ed/ii:', ed/ii)
                print('ed/ea:', ed/ea)
                print('Percentage Error:', np.absolute((Vp - Vpp)/Vp) * 100)
                Vp = (Vp + Vpp)/2
                ne = nee
                if np.absolute((Vp - Vpp)/Vp) > tol_iter:
                    continue
                else:
                    Vp = Vpp
                    ne = nee
                    print('')
                    print('__________End_Boltzmann_Solver__________')
                    print('')
                    print('Actual_Vp(eV):', Vpp)
                    print("Actual ne(1/m3): {:.2e}".format(nee))
                    Te = np.multiply(0.66666666666666666666666, np.trapz(((E ** 1.5)*z1), E))
                    np.savetxt('EEPF_4Pa_np5e15.txt', np.column_stack([E,z1]), delimiter = ',')
                    print('Actual_Te(eV): ', Te)       
                    print('ed/ii:', ed/ii)
                    print('ed/ea:', ed/ea)
                    plt.plot(E,z1)
                    fig = plt.gcf()
                    plt.grid(True)
                    fig.set_size_inches(10,10)
                    plt.yscale('log')
                    plt.ylim(10 ** -2.0, 10 ** -1)
                    plt.xlim(-1,30)
                    plt.xlabel('Energy (eV)', fontsize=15)
                    plt.ylabel('$EEPF^{-1.5}$ (eV)', fontsize=14)
                    
                    return (ne, Vp, Te, ed/ii, ed/ea, z1)
                break
            break
#Day_One = Boltzmann_Solver(EN_guess, Te_guess, ni_guess, 'Ar').EEPF_Vp_ne()
''' ########################### Star Dust #################################'''   
def ni_solver(Te_g, EN_g, bb_g, species_name):
    for nii in ni_range:
        ni = ni_guess + nii
        Box = Boltzmann_Solver(EN_g, Te_g, ni,bb_g, species_name).EEPF_Vp_ne()
        Power_Box = Power_Balance(Box[5], Box[0],ni, Box[1], Box[2], species_name).Net_Power()
        ParticleB = Particle_Balance(Box[5],species_name).Net_Production_Rate(Box[0],ni,Box[1],Box[2])
        print('')
        print('__________ni_Solver__________')
        print('')
        print('New ni(1/m3):', ni)
        print('Net Power:', Power_Box)
        print('Particle Balance:', ParticleB)
        if (abs(Power_Box) > tol_iter):
            continue 
        else:
            New_ni = ni_guess + nii
            New_Box = Boltzmann_Solver(EN_g, Te_g, New_ni, bb_g, species_name).EEPF_Vp_ne()
            New_Power_Box = Power_Balance(New_Box[5], New_Box[0], New_ni,
                                          New_Box[1], New_Box[2], species_name).Net_Power()
            New_ParticleB = Particle_Balance(New_Box[5],species_name).Net_Production_Rate(New_Box[0],
                                            New_ni,New_Box[1],New_Box[2])
            New_Particlei = Particle_Balance(New_Box[5],species_name).Ion_Generation(New_Box[0])
            Loss = Power_Balance(New_Box[5], New_Box[0], New_ni, New_Box[1], New_Box[2], species_name)
            print('')
            print('__________End_ni_Solver__________')
            print('')
            print('Electron Density (1/m3):', New_Box[0])
            print('Ion Density(1/m3):', New_ni)
            print('Particle Potential(eV):', New_Box[1]) 
            print('Electron Temperature(eV):', New_Box[2])
            print('Electron Temperature(eV):', New_Box[2])
            print('Momentum Transfer Power Loss(W):', Loss.L_momentum())
            print('Excitation Power Loss(W):', Loss.L_excitation())
            print('Ionization Power Loss(W):', Loss.L_ionization())
            print('Diffusion Power Loss(W):', Loss.L_ion_acc())
            print('Electron Collection Power Loss(W):', Loss.L_ecollection())
            print('Total Power Loss(W):', Loss.Power_Loss())
            print('Delta (ed/ii):', New_Box[3])
            print('Delta (ed/ea):', New_Box[4])
            print('Net Power:', New_Power_Box)
            print('Particle Balance:', New_ParticleB)
            print('Ion Gen Rate(1/s):', New_Particlei)
            return New_Box, New_ni
        break  
def EN_loop(Te_g, EN_g, bb_g, species_name):
    for ENN in EN_range:
        EN = EN_g + ENN
        Cube = ni_solver(Te_g, EN, bb_g, species_name)
        Particle_Cube = Particle_Balance(Cube[0][5], species_name).Net_Production_Rate(Cube[0][0],
                                        Cube[1],Cube[0][1],Cube[0][2])
        Power_Cube = Power_Balance(Cube[0][5],Cube[0][0],Cube[1],
                                   Cube[0][1], Cube[0][2], species_name).Net_Power()
        print('')
        print('__________EN_Solver__________')
        print('')
        print('New EN(Td):', EN)
        print('Power Balance:', Power_Cube)
        print('Particle Balance:', Particle_Cube)
        if (abs(Particle_Cube) > tol_iter):
            continue 
        else:
            New_EN = EN_g + ENN
            New_Cube = ni_solver(Te_g, New_EN, bb_g, species_name)
            New_Particle_Cube = Particle_Balance(New_Cube[0][5], species_name).Net_Production_Rate(New_Cube[0][0],
                                                New_Cube[1], New_Cube[0][1],New_Cube[0][2])
            New_Power_Cube = Power_Balance(New_Cube[0][5], New_Cube[0][0],
                                           New_Cube[1], New_Cube[0][1], New_Cube[0][2], species_name).Net_Power()
            LOS = Power_Balance(New_Cube[0][5], New_Cube[0][0], New_Cube[1],
                                New_Cube[0][1], New_Cube[0][2], species_name)
            New_C = abs((4 * np.pi * co.epsilon_0 * New_Cube[0][1] * rp) / co.e)
            print('')
            print('__________End_EN_Solver__________')
            print('')
            print('Power Balance:', New_Power_Cube)
            print('Particle Balance:', New_Particle_Cube)
            print('EN(Td):', New_EN)
            print('Electron Density(1/m3):', New_Cube[0][0])
            print('Ion Density(1/m3):', New_Cube[1])
            print('Particle Potential(eV):', New_Cube[0][1]) 
            print('Electron Temperature(eV):', New_Cube[0][2])
            print('Momentum Transfer Power Loss(W):', LOS.L_momentum())
            print('Excitation Power Loss(W):', LOS.L_excitation())
            print('Ionization Power Loss(W):', LOS.L_ionization())
            print('Ion Acceleration Power Loss(W):', LOS.L_ion_acc())
            print('Electron Collection Power Loss(W):', LOS.L_ecollection())
            print('Total Power Loss(W):', LOS.Power_Loss())
            print('Removal Energy(eV):', removal_energy(New_C))
            print('Hump Energy(eV):', removal_energy(New_C) + abs(New_Cube[0][1])) 
            print('Electron Temperature(eV):', New_Cube[0][2])
            print('Delta (ed/ii):', New_Cube[0][3])
            print('Delta (ed/ea):', New_Cube[0][4])
            print('Particle Charge:', New_C)
            plt.plot(E, New_Cube[0][5])
            fig = plt.gcf()
            plt.grid(True)
            fig.set_size_inches(10,10)
            plt.yscale('log')
            plt.ylim(10 ** -2, 10 ** -1)
            plt.xticks(np.arange(0, 16, step=1))
            plt.xlim(-0.2,16)
            plt.xlabel('Energy (eV)', fontsize=15)
            plt.ylabel('$EEPF^{-1.5}$ (eV)', fontsize=14)
            return New_Cube, New_EN
        break
#Day_One = EN_loop(Te_guess, EN_guess, bb_guess,  'Ar')
def Tp_loop(Te_g, EN_g, bb_g, species_name, species_one, species_two):
    for amp in amp_range:
        bb = bb_g + amp
        prism = EN_loop(Te_g, EN_g, bb_g, species_name)
        prism_C = abs((4 * np.pi * co.epsilon_0 * prism[0][0][1] * rp) / co.e)
        heat_box = Heat(prism_C, prism[0][0][2], prism[0][0][5], prism[0][1], prism[0][0][3], prism[0][0][4])
        temp_star = heat_box.Temp_Solver_Copper()[0]
        heat_prism = heat_box.thermion(temp_star)
        resonance_box = Resonance(prism[0][0][0], prism[0][1], prism[0][0][5], prism_C, species_one)
        photo_prism = resonance_box.photoion()
        quench_prism = resonance_box.quenchion()
        total_prism = heat_box.detachion()
        ragnarok_prism = heat_prism + photo_prism + quench_prism
        ghost_prism = abs((total_prism - ragnarok_prism)/total_prism)
        particle_prism = Particle_Balance(prism[0][0][5], species_name).Net_Production_Rate(prism[0][0][0],
                                         prism[0][1],prism[0][0][1],prism[0][0][2])
        power_prism = Power_Balance(prism[0][0][5],prism[0][0][0],prism[0][1], prism[0][0][1], 
                                    prism[0][0][2], species_name).Net_Power()
        print('')
        print('__________Tp_Solver__________')
        print('')
        print('Ghost Emission Ratio:',ghost_prism)
        print('therm/ii:', heat_prism)
        print('photoemission/ii:', photo_prism)
        print('quench/ii:', quench_prism)
        print('Total Loss/ ii', heat_prism + quench_prism + photo_prism)
        print('Total/ii', total_prism )
        print('Amplitude:', bb)
        print('New Tp:', temp_star)
        print('Power Balance:', power_prism)
        print('Particle Balance:', particle_prism)
        if ghost_prism > tol_iter:
            continue 
        else:
            new_bb = bb_g + amp
            new_prism = EN_loop(Te_g, EN_g, new_bb, species_name)
            new_prism_C = abs((4 * np.pi * co.epsilon_0 * new_prism[0][0][1] * rp) / co.e)
            new_heat_box = Heat(new_prism_C, new_prism[0][0][2], new_prism[0][0][5], 
                                new_prism[0][1], new_prism[0][0][3], new_prism[0][0][4])
            new_temp_star = new_heat_box.Temp_Solver_Copper()[0]
            new_heat_prism = new_heat_box.thermion(new_temp_star)
            new_resonance_box = Resonance(new_prism[0][0][0], new_prism[0][1], new_prism[0][0][5], new_prism_C, species_one)
            new_photo_prism = new_resonance_box.photoion()
            new_quench_prism = new_resonance_box.quenchion()
            new_total_prism = new_heat_box.detachion()
            new_ragnarok_prism = new_heat_prism + new_photo_prism + new_quench_prism
            new_ghost_prism = abs((new_total_prism - new_ragnarok_prism)/new_total_prism)
            new_particle_prism = Particle_Balance(new_prism[0][0][5], species_name).Net_Production_Rate(new_prism[0][0][0],
                                                 new_prism[0][1],new_prism[0][0][1],new_prism[0][0][2])
            new_power_prism = Power_Balance(new_prism[0][0][5],new_prism[0][0][0],
                                            new_prism[0][1], new_prism[0][0][1], new_prism[0][0][2], species_name).Net_Power()
            plos = Power_Balance(new_prism[0][0][5],new_prism[0][0][0],
                                 new_prism[0][1], new_prism[0][0][1], new_prism[0][0][2], species_name)
            L_relax = EERL(E, species_one, species_two, f_maxmin, new_prism[0][0][2],
                           new_prism[0][0][0], new_prism[0][1],new_prism[0][0][1])
            print('')
            print('__________End_Tp_Solver__________')
            print('')
            print('Ghost Emission Ratio:', new_ghost_prism)
            print('Power Balance:', new_power_prism)
            print('Particle Balance:', new_particle_prism)
            print('Particle Density(1/m3):', npp)
            print('Ion Density(1/m3):', new_prism[0][1])
            print('Electron Density(1/m3):', new_prism[0][0][0])
            print('Electron Temperature(eV):', new_prism[0][0][2])
            print('Momentum Transfer Power Loss(W):', plos.L_momentum())
            print('Excitation Power Loss(W):', plos.L_excitation())
            print('Ionization Power Loss(W):', plos.L_ionization())
            print('Ion Acceleration Power Loss(W):', plos.L_ion_acc())
            print('Electron Collection Power Loss(W):', plos.L_ecollection())
            print('EN(Td):', new_prism[1])
            print('Particle Charge:', new_prism_C)
            print('Particle Temperature(K):', new_temp_star)
            print('Hump Energy(eV):', (E_fermi + abs(new_prism[0][0][1])))
            print('Particle Potential(eV):', new_prism[0][0][1]) 
            print('Removal Energy(eV):', removal_energy(new_prism_C))
            print('Total Power Loss(W):', plos.Power_Loss())
            print('Total Power In(W):', P_in)
            print('Delta (ed/ii):', new_prism[0][0][3])
            print('Delta (ed/ea):', new_prism[0][0][4])
            print('therm/ii:', new_heat_prism)
            print('photoemission/ii:', new_photo_prism)
            print('quench/ii:', new_quench_prism)
            print('Amplitude of detach Cross Section:', new_bb)
            print('Pressure(Pa):', p)
            print('Sheath Voltage(V):', Vsh)
            plt.subplot(211)
            plt.plot(E, new_prism[0][0][5],linewidth=3)
            plt.plot(E, new_prism[0][0][5],linewidth=3)
            fig = plt.gcf()
            plt.grid(True)
            fig.set_size_inches(10,10)
            plt.yscale('log')
            plt.ylim(10 ** -2, 0.6*10 ** -1)
            plt.xticks(np.arange(0, 16, step=1))
            plt.xlim(-0.2,16)
            plt.xlabel('Energy (eV)', fontsize=15)
            plt.ylabel('$EEPF^{-1.5}$ (eV)', fontsize=14)
            plt.subplot(212)
            plt.plot(E, L_relax[0], linewidth=3, color= 'r')
            #plt.plot(E, L_relax[1], linewidth=3)
            plt.xlabel('Energy (eV)')
            plt.ylabel('Electron Energy Relaxation Length [m]')
            #plt.legend(['With Particles', 'Without Particles'])
            fig = plt.gcf()
            plt.grid(True)
            fig.set_size_inches(12,12)
            plt.rcParams.update({'font.size': 14})
            plt.xlim(0,5)
            pog = new_heat_box
            f_rad = pog.rad_copper(new_temp_star) / pog.Total_Loss_Copper(new_temp_star)
            f_cond = pog.conduction(new_temp_star) / pog.Total_Loss_Copper(new_temp_star)
            f_therm = pog.Thermionic() / pog.Total_Loss_Copper(new_temp_star)
            f_eke = pog.EKE() / pog.Total_Loss_Copper(new_temp_star)
            f_ike = pog.IKE()/ pog.Total_Loss_Copper(new_temp_star)
            f_e3p2 = pog.E_3p2()/ pog.Total_Loss_Copper(new_temp_star)
            f_e3p0 = pog.E_3p0()/ pog.Total_Loss_Copper(new_temp_star)
            f_ad = pog.AD()/ pog.Total_Loss_Copper(new_temp_star)
            f_atomh = pog.AtomicH()/ pog.Total_Loss_Copper(new_temp_star)
            f_las = pog.Laser_Heating()/ pog.Total_Loss_Copper(new_temp_star)
            f_others = f_ad + f_e3p2 + f_e3p0
            l_labels = 'Radiation', 'Conduction', 'Thermionic'
            l_sizes = [f_rad, f_cond, f_therm]
            print('------------------------------')
            print('total power loss (W)', pog.Total_Loss_Copper(new_temp_star))
            print('Radiation   %', f_rad*100)
            print('Conduction %', f_cond*100)
            print('Thermionic Emission %', f_therm*100)
            print('Electron Impaction %', f_eke*100)
            print('Ion Impaction %', f_ike*100)
            print('H3+ Asso-Disso %', f_ad*100)
            print('Atomic H2 Impaction %', f_atomh*100)
            print('Argon 3P2 Impaction %', f_e3p2*100)
            print('Argon 3P0 Impaction %',f_e3p0*100)
            print('------------------------------')
            h_labels = 'EKE', 'IKE', 'Atomic H', 'Others'
            h_sizes = [f_eke, f_ike, f_atomh, f_others]
            l_explode = (0,0,0)
            h_explode = (0,0,0.0,0)
            cmap = plt.cm.prism
            fig = plt.figure(figsize=(20,20))
            cmap = plt.cm.prism                    
            fig, (ax1, ax2) = plt.subplots(1,2)
            fig.suptitle('Particle Temperature: {} K @ {} Pa @ np: 1E15 (1/m3)'.format(round(new_temp_star),p), fontweight='bold', y=0.9)
            plt.subplots_adjust(wspace= 100, hspace= 100, top=0.7, bottom=0, left=0.1, right=0.9)
            ax1.pie(h_sizes, explode=h_explode, autopct='%1.1f%%',
                    radius=85, shadow=True, startangle=90, labels = h_labels,
                    labeldistance=1.1)
            ax2.pie(l_sizes, explode=l_explode,labels= l_labels, autopct='%1.1f%%', 
                    radius=85, shadow=True, startangle=90, labeldistance=1.1)
            '''
            temperature = np.linspace(300,2600,10000)
            Con = new_heat_box.conduction(temperature)
            Radiation = new_heat_box.rad_copper(temperature)
            Term = new_heat_box.Thermionic()*(temperature/temperature)
            ionkin = new_heat_box.IKE()*(temperature/temperature)
            elekin = new_heat_box.EKE() * (temperature/temperature)
            e3p2kin = new_heat_box.E_3p2()*(temperature/temperature)
            e3p0kin = new_heat_box.E_3p0()*(temperature/temperature)
            adion = new_heat_box.AD()*(temperature/temperature)
            atomic = new_heat_box.AtomicH()*(temperature/temperature)
            lash = new_heat_box.Laser_Heating()*(temperature/temperature)
            plt.subplot(313)
            plt.plot(temperature, Con, linewidth=3)
            plt.plot(temperature, Radiation, linewidth=3)
            plt.plot(temperature, Term, linewidth=3)
            plt.plot(temperature, ionkin, linewidth=3)
            plt.plot(temperature, elekin, linewidth=3)
            plt.plot(temperature,e3p2kin, linewidth=3 )
            plt.plot(temperature,e3p0kin, linewidth=3)
            plt.plot(temperature, adion, linewidth=3)
            plt.plot(temperature, atomic, linewidth=3)
            plt.plot(temperature, lash, linewidth=3)
            plt.axvline(x=new_temp_star,color='r', linestyle='--')
            plt.xlabel('Temperature [K]')
            plt.xlim(300,2600)
            plt.ylabel('Heat Loss [W]')
            plt.legend(['Conduction','Radiation', 'Thrmionic Emission', 'IKE',
                        'EKE', 'E_3p2','E_3P0', 'Asso/Disso:ciation H3+',
                        'Radical atomic H', 'Laser Heating'])
            '''
            #np.savetxt('EEPF_4Pa_dist.txt',(new_prism[0][0][5]), fmt='%.15f')
            #np.savetxt('EERL_Wo_2.5pa.txt',L_relax[1], fmt='%.15f')
            #np.savetxt('EEPF_WF4.7_2.5pa.txt',L_relax[1], fmt='%.15f')
            return new_prism, new_temp_star, new_bb
        break
Day_One = Tp_loop(Te_guess, EN_guess, bb_guess,  'H2', 'Ar', 'H2')
            
            










