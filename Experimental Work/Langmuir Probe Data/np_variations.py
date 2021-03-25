import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as co
#path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/01-31-Cu_double_Ar/Ar Addition/ne_ni_Te_hu.txt'
#path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/12-19-copper/pressure study/ne_ni_Te_hu.txt'
#T = pd.read_csv(path_dir, delimiter= ',')
#npp = np.array(T.iloc[0:8,0]).astype(float)
ni = np.array([5.4182850127153837e+17,5.575845067809493e+17,5.644985321522464e+17,
               5.633676962890809e+17,5.59653141903224e+17,5.579332919403409e+17])
ne = np.array([1.279906940691818e+16,1.3606553138772766e+16,1.3886957663548376e+16,
               1.408202574044621e+16,1.4046611190245878e+16,1.3949349097469878e+16])
Te = np.array([4.97463708631586,5.098855971950952,5.273347310076409,
               5.470563068244049,5.706148238154203,5.966881334713291])
Pin = np.array([46.5, 56.7, 66.1, 75.5, 84.6, 95.5])
#print(ni)
rp = 16e-9

c_lim = (4*np.pi* co.epsilon_0*4.7/co.e)*rp + (3/8)
print('Charge limit', c_lim)

'''
lmom = np.array(T.iloc[0:8,4]).astype(float)
lex = np.array(T.iloc[0:8,5]).astype(float)
lion = np.array(T.iloc[0:8,6]).astype(float)
lacc = np.array(T.iloc[0:8,7]).astype(float)
lecoll = np.array(T.iloc[0:8,8]).astype(float)
EN = np.array(T.iloc[0:8,9]).astype(float)
ch = np.array(T.iloc[0:8,10]).astype(float)
Tp = np.array(T.iloc[0:8,11]).astype(float)
hump = np.array(T.iloc[0:8,12]).astype(float)
Vp= np.array(T.iloc[0:8,13]).astype(float)
Erem = np.array(T.iloc[0:8,14]).astype(float)
Pin = np.array(T.iloc[0:8,16]).astype(float)
edii = np.array(T.iloc[0:8,17]).astype(float)
edea = np.array(T.iloc[0:8,18]).astype(float)
bb = np.array(T.iloc[0:8,19]).astype(float)
p = np.array(T.iloc[0:8,20]).astype(float)
plt.title('Power Study @ W: 4.7 eV, R: 30 nm, H Heating')
wf = np.array(T.iloc[0:8,22]).astype(float)
'''


'''
plt.plot(Pin,lmom, linewidth=4)
plt.plot(Pin,lex, linewidth=4)
plt.plot(Pin,lion, linewidth=4)
plt.plot(Pin,lacc, linewidth=4)
plt.plot(Pin,lecoll, linewidth=4)
#plt.plot(Pin,Pin)
fig = plt.gcf()
plt.grid(True)
fig.set_size_inches(10, 10)
plt.xlabel('Power In (W)', fontsize=14)
plt.ylabel('Power Loss (W)', fontsize=14)
plt.rcParams.update({'font.size': 9})
plt.legend(['L_momentum', 'L_excitation', 'L_ionization', 'L_ion_accel', 'L_e_collec (Pin: 90W)'])
#plt.ylim(-1,80)
'''




def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


fig, host = plt.subplots()
fig.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()

# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
par2.spines["right"].set_position(("axes", 1.2))
# Having been created by twinx, par2 has its frame off, so the line of its
# detached spine is invisible.  First, activate the frame but make the patch
# and spines invisible.
make_patch_spines_invisible(par2)
# Second, show the right spine.
par2.spines["right"].set_visible(True)

'''
p1, = host.plot(Pin,Te, "b-", label="Electron Temperature (eV)", linewidth=4)
p2, = par1.plot(Pin,EN, "r-", label="EN (Td)", linewidth=4)
p3, = par2.plot(Pin, ch, "g-", label="Particle Charge", linewidth=4)
host.set_xlabel("Power (W) ")
host.set_ylabel("Electron Temperature (eV)")
par1.set_ylabel("EN (Td)")
par2.set_ylabel("Particle Charge")
'''

#fig = plt.gcf()
#plt.grid(True)
fig.set_size_inches(12, 12)
#plt.title('100 ks')

#plt.xscale('log')
#plt.ylim(10 ** -3.0, 10 ** 0)
#plt.xlim(-1,13)
plt.rcParams.update({'font.size': 20})

host.set_xlim(35, 72)
host.set_ylim(1e9,5e11)
host.set_yscale('log')
par1.set_ylim(1e9,5e11)
par1.set_yscale('log')

par2.set_ylim(5.2,7.2)

'''
p1, = host.plot(Pin,Erem, "m-", label="Removal Energy (eV)", linewidth=4)
p2, = par1.plot(Pin,-Vp, "k-", label="|Particle Potential| (eV)", linewidth=4)
p3, = par2.plot(Pin, hump, "--", label="Peak Energy (eV)", linewidth=4)
host.set_xlabel("Power (W)")
host.set_ylabel("Removal Energy (eV)")
par1.set_ylabel("|Particle Potential| (eV)")
par2.set_ylabel("Peak Energy (eV)")
'''

p1, = host.plot(Pin,ne*1e-6*4.070804524135282, color="b", label="$n_e({cm^{-3}}$)", linewidth=4)
p2, = par1.plot(Pin,ni*1e-6, color="g", label="$n_i({cm^{-3}}$)", linewidth=4)
p3, = par2.plot(Pin,Te,color="r", label="$T_e(K)$", linewidth=4)
host.set_xlabel("Pressure (mTorr)")
host.set_ylabel("$n_e({cm^{-3}}$)")
par1.set_ylabel("$n_i({cm^{-3}}$)")
par2.set_ylabel("$T_e(K)$")

print((ni -(ne*4.070804524135282))/c_lim)

'''
p1, = host.plot(Pin,edii, "b-", label=" ed/ ii", linewidth=4)
p2, = par1.plot(Pin,edea, "g-", label="ed/ ea", linewidth=4)
p3, = par2.plot(Pin,bb, "r-", label="Amplitude of Detachment CS", linewidth=4)
host.set_xlabel("Power (W)")
host.set_ylabel("Electron Detach Freq / Ion Attach Freq")
par1.set_ylabel("Electron Detach Freq / Electron Attach Freq)")
par2.set_ylabel("Amplitude of Detachment CS")
'''



host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)

host.tick_params(axis='y', direction= 'in',colors=p1.get_color(), **tkw)
par1.tick_params(axis='y', direction= 'in', colors=p2.get_color(), **tkw)
par2.tick_params(axis='y', direction= 'in', colors=p3.get_color(), **tkw)
host.tick_params(axis='x', direction= 'in' ,**tkw)

lines = [p1, p2, p3]

host.legend(lines, [l.get_label() for l in lines])
plt.savefig('pressure_final.png', dpi = 600)
