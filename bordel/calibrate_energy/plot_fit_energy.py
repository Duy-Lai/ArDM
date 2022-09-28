import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

file = '/home/duy/polybox/master/ardm/calibrate_energy/fit_data.csv'
data = pd.read_csv(file, sep=',')

# Plotting parameters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size + 2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams["legend.markerscale"] = 1
plt.rcParams['axes.axisbelow'] = True
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)

energy = 0.5*(data['start']+data['end'])
pts = np.linspace(0,100,100)

def linear_fit(x, a, b):
        return a*x+b
        
popt, _ = curve_fit(linear_fit, energy[2:], data['mean_F90'][2:])

fig, ax = plt.subplots()
p1 = plt.errorbar(energy, data['mean_F90'], yerr=data['sigma_F90'],
                xerr=0.5*(data['end']-data['start']), fmt='o', markersize=6, capsize=4)
p2, = ax.plot(pts, linear_fit(pts, *popt), 'r--', zorder=0)
# plt.plot(0.5*(data['start']+data['end']), data['mean_F90'], 'x-')
# plt.fill_between(0.5*(data['start']+data['end']),
        # data['mean_F90']-data['sigma_F90'], data['mean_F90']+data['sigma_F90'], alpha=0.25)
lgnd = ax.legend([p1,p2], [r'2D Gaussian fit results',r'Linear fit'], framealpha=0.5)
plt.grid()
plt.xlim(0.,100.)
plt.ylim(0.,1.)
plt.xlabel('S1 [p.e.]')
plt.ylabel('F90')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.grid(alpha=0.5, zorder=0)
fig.tight_layout()
plt.savefig('/home/duy/polybox/master/ardm/calibrate_energy/fit_mean_F90_vs_S1L.png', dpi=300)

popt, _ = curve_fit(linear_fit, energy, data['mean_logS2maxLoS1L'])

fig, ax = plt.subplots()
p1 = ax.errorbar(0.5*(data['start']+data['end']), data['mean_logS2maxLoS1L'], yerr=data['sigma_logS2maxLoS1L'],
                xerr=0.5*(data['end']-data['start']), fmt='o', markersize=6, capsize=4)
p2, = ax.plot(pts, linear_fit(pts, *popt), 'r--', zorder=0)
lgnd = ax.legend([p1,p2], [r'2D Gaussian fit results',r'Linear fit'], framealpha=0.5)
plt.grid()
plt.xlim(0.,100.)
plt.ylim(-3.,5)
plt.xlabel('S1 [p.e.]')
plt.ylabel(r'$\log\left(\frac{S2}{S1}\right)$')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.grid(alpha=0.5, zorder=0)
fig.tight_layout()
plt.savefig('/home/duy/polybox/master/ardm/calibrate_energy/fit_mean_logS2S1_vs_S1L.png', dpi=300)
