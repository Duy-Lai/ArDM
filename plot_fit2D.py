# scripts similar to plot_warios.py, to make all the plots in the .root file produced by Fit2D.cpp
# also requires the CM fit results, CM_results.csv, produced by Fit2D.cpp
# written by Duy Lai, 09.2022

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors, patches
from matplotlib.ticker import (MultipleLocator, LogLocator, ScalarFormatter)
import uproot, os

slice = 5

root_file = 'fit2D/histos.root'
histo = np.array(uproot.open(root_file).keys())

tab = pd.read_csv('neutron/CM_results.csv')

def format_number(num):
    x = int(np.log10(abs(num)))
    if num < 1:
        x -= 1
    if x <= -3:
        return r'{:>.2e}'.format(num)
    else:
        return (r'{:>.' + str(3-x) + r'f}').format(num)

def StatsBox():
    x = 0.5 * (data[1][:-1] + data[1][1:])
    y = 0.5 * (data[2][:-1] + data[2][1:])
    y, x = np.meshgrid(y, x)
    mean = np.array([np.average(x, weights=data[0]), np.average(y, weights=data[0])])
    std = np.sqrt(np.array([np.average((x-mean[0])**2, weights=data[0]), np.average((y-mean[1])**2, weights=data[0])]))
    handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 5
    labels = []
    if tcontent[0][:-1] == 'CM_Neutron_ExtPdf':
        labels.append(r'Rec. ER events    {:,.0f}'.format(tab.N_ER_Neutron[slice]))
        labels.append(r'Rec. NR events    {:,.0f}'.format(tab.N_NR_Neutron[slice]))
    elif tcontent[0][:-1] == 'CM_DP_ExtPdf':
        labels.append(r'Rec. ER events    {:,.0f}'.format(tab.N_ER_DP[slice]))
        labels.append(r'Rec. NR events    {:,.0f}'.format(tab.N_NR_DP[slice]))
    elif tcontent[0][:-1] in ['histo_DP','histo_Neutron']:
        labels.append(r'# of events    {:,.0f}'.format(data[0].sum()))
    else:
        labels.append(r'Rec. # of events    {:,.0f}'.format(data[0].sum()))
        # labels.append(r'Mean x           ' + format_number(mean[0]))
        # labels.append(r'Mean y           ' + format_number(mean[1]))
        # labels.append(r'Std Dev x       ' + format_number(std[0]))
        # labels.append(r'Std Dev y       ' + format_number(std[1]))
    ax.legend(handles, labels, loc='best', fontsize= font_size+mag-2, 
          fancybox=True, framealpha=0.5, 
          handlelength=0, handletextpad=0)

mag = 6

# Plotting parameters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size + 2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size - 4
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams["legend.markerscale"] = 3.5
plt.rc('xtick', labelsize=font_size+mag)
plt.rc('ytick', labelsize=font_size+mag)

# for n in [149]:
for n in np.arange(histo.shape[0]):
    print('{} -- {}'.format(n+1, histo[n][:-2]))
    tcontent = []
    k = 0
    for j in np.arange(len(histo[n])):
        if histo[n][j] == ' ' or j == len(histo[n])-1:
            tcontent.append(histo[n][k:j])
            k = j + 1
    if 1 + 1 == 2:
    # if float(tcontent[2][1:-1]) == slice * 10:
        fig, ax = plt.subplots(figsize=(7,5))
        data = uproot.open(root_file)[histo[n]].to_numpy()
        Y, X = np.meshgrid(data[2], data[1])
        Z = data[0]
        low_limit = 8 / np.power(data[1].shape[0]-1, 2)
        for i in np.arange(Z.shape[0]):
            for j in np.arange(Z.shape[1]):
                if tcontent[0][:-1] in ['Neutron_discriminator', 'Electron_discriminator', 'CM_Neutron_error', 'CM_DP_error']:
                    if Z[i][j] == 0.0:
                        Z[i][j] = np.NaN
                else:
                    if Z[i][j] < low_limit:
                        Z[i][j] = 0.0
        if tcontent[0][:-1] in ['Neutron_discriminator', 'Electron_discriminator', 'CM_Neutron_error', 'CM_DP_error']:
            plt.pcolor(X, Y, Z, vmin=-1.0, vmax=1.0, cmap='jet')
            plt.colorbar()
        else:
            plt.pcolor(X, Y, Z, norm=colors.LogNorm(vmin=low_limit), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            # StatsBox()
        plt.title(r'{:.0f}$\leq$S1LY<{:.0f} p.e.'.format(float(tcontent[2][1:-1]), float(tcontent[3][:-2])), fontsize=font_size+mag+6)
        # plt.yticks(np.arange(-3,6))
        ax.grid(alpha=0.5, zorder=0)
        # plt.title(histo[i][:-2])
        plt.xlabel(r'S1F90', fontsize=font_size+mag+8)
        plt.ylabel(r'S2maxLos1Lcorr', fontsize=font_size+mag+8)
        fig.tight_layout()
        try:
            plt.savefig('thesis_plots/fit2D/{}/{:.0f}_{:.0f}.png'.format(tcontent[0][:-1], float(tcontent[2][1:-1]), float(tcontent[3][:-2])), dpi=300)
        except FileNotFoundError:
            try:
                os.mkdir('thesis_plots/fit2D/{}'.format(tcontent[0][:-1]))
            except FileNotFoundError:
                os.mkdir('thesis_plots/fit2D')
                os.mkdir('thesis_plots/fit2D/{}'.format(tcontent[0][:-1]))
            plt.savefig('thesis_plots/fit2D/{}/{:.0f}_{:.0f}.png'.format(tcontent[0][:-1], float(tcontent[2][1:-1]), float(tcontent[3][:-2])), dpi=300)
        # plt.savefig('thesis_plots/fit2D/{}'.format(tcontent[0][:-1]), dpi=300)
        plt.close()