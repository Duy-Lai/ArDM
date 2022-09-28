import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
from scipy.optimize import curve_fit
import uproot

file = '/home/duylai/ardm/calibrate_energy_sniper/fit_data.csv'
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
        
popt, _ = curve_fit(linear_fit, energy[1:], data['mean_F90'][1:])

fig, ax = plt.subplots()
histo = uproot.open('/mnt/raid/users/duylai/warios_sniper_NR/SuperWario.root')['pre_F90_vs_S1L']
histo_data = histo.to_numpy()
X, Y = np.meshgrid(histo_data[2], histo_data[1])
Z = histo_data[0]
ax.pcolor(Y, X, Z, norm=colors.LogNorm(), cmap='jet', zorder=0)
p1 = plt.errorbar(energy, data['mean_F90'], yerr=data['sigma_F90'],
                xerr=0.5*(data['end']-data['start']), c='k', fmt='o', markersize=6, capsize=4, zorder=1)
p2, = ax.plot(pts, linear_fit(pts, *popt), 'r--')
# plt.plot(0.5*(data['start']+data['end']), data['mean_F90'], 'x-')
# plt.fill_between(0.5*(data['start']+data['end']),
        # data['mean_F90']-data['sigma_F90'], data['mean_F90']+data['sigma_F90'], alpha=0.25)
lgnd = ax.legend([p1,p2], [r'2D Gaussian fit results',r'Linear fit (1st data point excluded)'], framealpha=0.5, loc='lower left')
plt.grid()
plt.xlim(0.,100.)
plt.ylim(0.,1.)
plt.xlabel('S1 [p.e.]')
plt.ylabel('F90')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.grid(alpha=0.5, zorder=0)
fig.tight_layout()
plt.savefig('/home/duylai/ardm/calibrate_energy_sniper/fit_mean_F90_vs_S1L.png', dpi=300)

popt, _ = curve_fit(linear_fit, energy, data['mean_logS2maxLoS1L'])
popt_up, _ = curve_fit(linear_fit, energy, data['mean_logS2maxLoS1L']+data['sigma_logS2maxLoS1L'])
popt_down, _ = curve_fit(linear_fit, energy, data['mean_logS2maxLoS1L']-data['sigma_logS2maxLoS1L'])

fig, ax = plt.subplots()
p1 = ax.errorbar(energy, data['mean_logS2maxLoS1L'], yerr=data['sigma_logS2maxLoS1L'],
                xerr=0.5*(data['end']-data['start']), fmt='o', markersize=6, capsize=4)
p2, = ax.plot(pts, linear_fit(pts, *popt), 'r--', zorder=0)
ax.plot(pts, linear_fit(pts, *popt_up), 'r:', zorder=0)
ax.plot(pts, linear_fit(pts, *popt_down), 'r:', zorder=0)
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
plt.savefig('//home/duylai/ardm/calibrate_energy_sniper/fit_mean_logS2S1_vs_S1L.png', dpi=300)
