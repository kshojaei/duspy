import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import scipy.stats as sps
from scipy.optimize import curve_fit
import glob
import seaborn as sns
from scipy.stats import norm



#########
#path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/02-27-cooper-H2/H2 Content/'
#path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/01-31-Cu_double_Ar/Ar Addition/'
path_dir = r'/Users/scienceman/Desktop/Project_Titan/Titan/02-27-cooper-H2/Pressure/'


T = pd.read_csv(path_dir+'10h2__46.5mTorr_EEPF.txt', delimiter= ',')
T2 = pd.read_csv(path_dir+'10h2__56.7mTorr_EEPF.txt', delimiter= ',')
T3 = pd.read_csv(path_dir+'10h2__66.1mTorr_EEPF.txt', delimiter= ',')
T4 = pd.read_csv(path_dir+'10h2__75.5mTorr_EEPF.txt', delimiter= ',')
T5 = pd.read_csv(path_dir+'10h2__84.6mTorr_EEPF.txt', delimiter= ',')
T6 = pd.read_csv(path_dir+'10h2__95.5mTorr_EEPF.txt', delimiter= ',')
#T7 = pd.read_csv(path_dir+'30Ar.txt', delimiter= ',')
#T8 = pd.read_csv(path_dir+'35Ar.txt', delimiter= ',')
#T9 = pd.read_csv(path_dir+'45Ar.txt', delimiter= ',')
#T10 = pd.read_csv(path_dir+'50Ar.txt', delimiter= ',')


E = np.array(T.iloc[:,0]).astype(float)
E2 = np.array(T2.iloc[:,0]).astype(float)
E3 = np.array(T3.iloc[:,0]).astype(float)
E4 = np.array(T4.iloc[:,0]).astype(float)
E5 = np.array(T5.iloc[:,0]).astype(float)
E6 = np.array(T6.iloc[:,0]).astype(float)
#E7 = np.array(T7.iloc[:,0]).astype(float)
#E8 = np.array(T8.iloc[:,0]).astype(float)
#E9 = np.array(T9.iloc[:,0]).astype(float)
#E10 = np.array(T10.iloc[:,0]).astype(float)






W_15 = np.array(T.iloc[:,1]).astype(float)
W_27 = np.array(T2.iloc[:,1]).astype(float)
W_37 = np.array(T3.iloc[:,1]).astype(float)
W_51 = np.array(T4.iloc[:,1]).astype(float)
W_65 = np.array(T5.iloc[:,1]).astype(float)
W_84 = np.array(T6.iloc[:,1]).astype(float)
#W_1 = np.array(T7.iloc[:,1]).astype(float)
#W_2 = np.array(T8.iloc[:,1]).astype(float)
#W_3 = np.array(T9.iloc[:,1]).astype(float)
#W_4 = np.array(T10.iloc[:,1]).astype(float)





fig, ax = plt.subplots()
fig = plt.gcf()
#plt.plot(E,W_15*1E-6, linewidth=4)
ax.plot(E2,W_27*1E-6, linewidth=4, color= 'navy')
#ax.plot(E3,W_37*1E-6, linewidth=4)
plt.plot(E4,W_51*1E-6, linewidth=4, color = 'gold')
plt.plot(E5,W_65*1E-6, linewidth=4, color= 'r')
plt.plot(E6,W_84*1E-6, linewidth=4, color= 'lime')
#plt.plot(E7,W_1*1E-6, linewidth=4)
#plt.plot(E8,W_2*1E-6, linewidth=4)
#plt.plot(E9,W_3*1E-6, linewidth=4)
#plt.plot(E10,W_4*1E-6, linewidth=4)
#fig = plt.gcf()
plt.grid(True)
ax.tick_params(direction= 'in')
plt.xlim(-0.5,16)

plt.text(0.1,0.8e8, r'(a)', fontsize=13)
plt.xlim(-0.5,16)


Pin = np.array([46.5, 56.7, 66.1, 75.5, 84.6, 95.5])

plt.tick_params(direction= 'in')
plt.tick_params(axis="y",direction="in")
plt.rcParams.update({'font.size': 13})
plt.ylim(1E7, 1.1E8)
plt.xlim(-0.5,16)
plt.xticks(np.arange(0, 18, step=2))
#plt.yticks(np.arange(1e8, 1e9, step=4))
#ax.text(2, 27, r'(a)', fontsize=13)
#ax.text(2, 27, r'(a)', fontsize=13)

plt.xlabel('Energy (eV)', fontsize=13)
plt.ylabel('$f_p(cm^{-3}eV^{-3/2})$', fontsize=13)
plt.grid(True)
#fig.set_size_inches(10, 10)
plt.legend(['56.7 mTorr', '75.5 mTorr', '84.6 mTorr', '95.5 mTorr'])
#ax.tick_params(direction= 'in')
#ax.tick_params(axis="y",direction="in")
plt.yscale('log')
plt.tick_params(direction='in')
plt.tick_params(which='minor', direction='in')


plt.savefig('pressure_eedf.png', dpi = 600)




'''

size = np.array(T.iloc[:,13]).astype(float)


kwargs = dict(hist_kws={'alpha':0.9})



fig, ax = plt.subplots()
sns.distplot(size, color="blue", kde=False, **kwargs)

plt.xlabel('Particle Size (nm)', fontsize=20, fontweight='bold')
plt.ylabel('Particle Fraction (%)', fontsize=20, fontweight='bold')
fig = plt.gcf()
plt.grid(True)
fig.set_size_inches(12,12)
plt.rcParams.update({'font.size': 20})

#sns.distplot(size, color="blue", hist=True)
#sns.kdeplot(size, bw=9, linewidth=4, color='red')


mu, sigma = norm.fit(size)
print(mu, sigma)
x = np.linspace(0,100,1000)
best_fit_line = norm.pdf(x, mu, sigma)*600
plt.plot(x, best_fit_line, color='red', linewidth=5)
plt.xticks(np.arange(0, 100, step=10))
plt.yticks(np.arange(0, 26, step=4))
ax.tick_params(direction= 'in')
ax.text(60, 3, r'AVG: 32.4 nm', fontsize=20, fontweight='bold')
ax.text(60, 2, r'SD: 11.3 nm ', fontsize=20, fontweight='bold')


plt.legend(['Gaussian Fit','Particle Size Distribution'])
plt.xlim(0,91)
plt.savefig('plot_name.png', dpi = 600)
# add a 'best fit' line



#density = np.sum(size, axis=0)
#density /= np.trapz(density, size)


#sns.distplot(size, color="red", label="Compact", **kwargs)

#sns.distplot(x2, color="orange", label="SUV", **kwargs)
#sns.distplot(x3, color="deeppink", label="minivan", **kwargs)

'''































