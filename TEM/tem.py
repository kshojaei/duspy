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
path_dir = r'/Users/scienceman/Desktop/Langmuir/TEM_COPPER_GOLD/c5/Copper_Plasma.csv'


T = pd.read_csv(path_dir, delimiter= ',')

size = np.array(T.iloc[:,13]).astype(float)


kwargs = dict(hist_kws={'alpha':0.9})



fig, ax = plt.subplots()
sns.distplot(size, color="blue", kde=False, **kwargs)
plt.rcParams.update({'font.size': 13})
plt.xlabel('Particle Size (nm)', fontsize=13)
plt.ylabel('Particle Fraction (%)', fontsize=13)
fig = plt.gcf()
plt.grid(True)
#fig.set_size_inches(12,12)


#sns.distplot(size, color="blue", hist=True)
#sns.kdeplot(size, bw=9, linewidth=4, color='red')


mu, sigma = norm.fit(size)
print(mu, sigma)
x = np.linspace(0,100,1000)
plt.rcParams.update({'font.size' : 13})  
best_fit_line = norm.pdf(x, mu, sigma)*550
gaussian = (1/(np.sqrt(2*np.pi)*sigma)) * np.exp(-((x- mu)**2)/(2*(sigma**2)))*1e-9*320
plt.plot(x,best_fit_line, color='red', linewidth=5)
plt.xticks(np.arange(0, 100, step=10))
plt.yticks(np.arange(0, 26, step=4))
ax.tick_params(direction= 'in')
#ax.text(2, 27, r'(a)', fontsize=13)
ax.text(63, 5, r'AVG: 32.4 nm', fontsize=13)
ax.text(63, 2, r'SD: 11.3 nm ', fontsize=13)
#ax.text(34, 21, r'Copper Nanoparticles ', fontsize=13, fontweight= 'bold')
plt.legend(['Gaussian Fit','Particle Size Distribution'])
plt.xlim(0,90)
plt.ylim(0,30)
#plt.rcParams["font.weight"] = "bold"
#plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams.update({'font.size': 13})
plt.savefig('plot_feb.png', dpi = 600)
# add a 'best fit' line



#density = np.sum(size, axis=0)
#density /= np.trapz(density, size)


#sns.distplot(size, color="red", label="Compact", **kwargs)

#sns.distplot(x2, color="orange", label="SUV", **kwargs)
#sns.distplot(x3, color="deeppink", label="minivan", **kwargs)








































