#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import (MultipleLocator, LogLocator, ScalarFormatter)
import uproot

# type = 'dp_ER'
type = 'neutron'
# type = 'sniper_NR'
# histo_name ='S2maxDT'
# histo_name = 'S1LY'
# histo_name = 'RoI_F90_vs_S1L'
histo_name = 'RoI_logS2maxLoS1L_vs_F90'
# histo_name = 'S1Nx_vs_nS2'
save_fig = 0

# if type[:2] == 'dp' and histo_name[:3] == 'pre':
    # type += '_small'

if save_fig == 0:
    show_fig = 1
else:
    show_fig = 0

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
plt.rcParams['axes.axisbelow'] = False
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)

def F90_lower(x):
    return 0.288+0.303*np.exp(-0.073*x)+5.0*(0.031+1.687/(x+3.805))

def F90_upper(x):
    return 0.737+0.101*np.exp(-0.026*x)+3.0*(0.238+345.822/(x-2139.788))

def RoI_cuts():
    xps = np.linspace(10,100,200)
    lower_cut = F90_lower(xps)
    upper_cut = F90_upper(xps)
    plt.plot(xps, lower_cut, 'k-.', linewidth=2)
    plt.plot(xps, upper_cut, 'k-.', linewidth=2)
    plt.plot([100,100], [lower_cut[-1],upper_cut[-1]], 'k-.', linewidth=2)
    plt.xlim(0,125)
    # plt.fill_between(xps, lower_cut, upper_cut, color='k', alpha=0.4)

if histo_name == 'S2maxDT':
    # recoil_type = 'NR'
    recoil_type = 'ER'
    scattering = ''
    types = ['dp_' + recoil_type + scattering + '/1bump/', 'neutron_' + recoil_type + scattering + '/1bump/',
                'dp_' + recoil_type + scattering + '/2bump/', 'neutron_' + recoil_type + scattering + '/2bump/']
    if recoil_type == 'ER':
        other_types = ['dp_NR' + scattering + '/1bump/', 'neutron_NR' + scattering + '/1bump/',
                'dp_NR' + scattering + '/2bump/', 'neutron_NR' + scattering + '/2bump/']
    elif recoil_type == 'NR':
        other_types = ['dp_ER' + scattering + '/1bump/', 'neutron_ER' + scattering + '/1bump/',
                'dp_ER' + scattering + '/2bump/', 'neutron_ER' + scattering + '/2bump/']
    fig, ax = plt.subplots(figsize=(7,5))
    hist_range = (0.0,1.4)
    density = False
    for i in np.arange(len(types)):
        root_file = '/mnt/raid/users/duylai/warios_' + types[i] + 'SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data = np.array(tree[histo_name], dtype=np.float64)
        if density == False:
            other_root_file = '/mnt/raid/users/duylai/warios_' + other_types[i] + 'SuperWario.root'
            other_tree = uproot.open(other_root_file + ':HebingTree')
            other_data = np.array(other_tree[histo_name], dtype=np.float64)
            weights = np.ones(data.shape[0])/(data.shape[0]+other_data.shape[0])
            # pre_histo = uproot.open(root_file)['pre_F90_vs_S1L']
            # sum_pre_events = 0
            # pre_data = pre_histo.to_numpy()
            # for i in np.arange(len(pre_data[1])-1):
                # for j in np.arange(len(pre_data[2])-1):
                    # if pre_data[1][i] <= 100:
                        # sum_pre_events += pre_data[0][i,j]
            # weights = np.ones(data.shape[0])/pre_histo.member('fEntries')
        else:
            weights = np.ones(data.shape[0])
        histo = plt.hist(data, range=hist_range, bins=200, weights=weights, log=0, density=density, histtype='step', label=type)
        S2maxFracS2L = np.array(tree['S2maxFracS2L'], dtype=np.float64)
    ax.grid(alpha=0.5, zorder=0)
    plt.xlim(hist_range)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    plt.xlabel(r'Drift time [ms]', fontsize=font_size+2)
    plt.ylabel(r'Fraction of total events (per 7 $\mu$s)')
    if recoil_type == 'NR':
        plt.ylim(0,5e-4)
        # plt.ylim(0,2e-5)
        plt.legend([r'$\log(S2/S1)\geq0.2$, DP', r'$\log(S2/S1)\geq0.2$, Neutron',
                r'$\log(S2/S1)<0.2$, DP', r'$\log(S2/S1)<0.2$, Neutron'], loc=0, fontsize=font_size)
        title_text = r'RoI'
    elif recoil_type == 'ER':
        plt.ylim(0,2e-2)
        # plt.ylim(0,1.2e-3)
        # ax.yaxis.set_major_locator(MultipleLocator(1e-6))
        plt.legend([r'$\log(S2/S1)\geq0.5$, DP', r'$\log(S2/S1)\geq0.5$, Neutron',
                r'$\log(S2/S1)<0.5$, DP', r'$\log(S2/S1)<0.5$, Neutron'], loc=0, fontsize=font_size)
        title_text = r'Inverse RoI'
    if density == True:
        plt.ylim(0,7)
        plt.ylabel(r'Density')
    if scattering == '_single':
        title_text += r', single scatters'
    elif scattering == '_multi':
        title_text += r', multiple scatters'
    plt.title(title_text, fontsize=font_size, color='r')
    fig.tight_layout()
    if save_fig:
        name = 'python/' + histo_name
        if density == True:
            name += '_density'
        if recoil_type == 'NR':
            plt.savefig(name + '_RoI' + scattering + '.png', dpi=300)
        elif recoil_type == 'ER':
            plt.savefig(name +'_inverse_RoI' + scattering + '.png', dpi=300)
    if show_fig:
        plt.show()

elif histo_name == 'S1LY':
    recoil_type = 'ER'
    fig, ax = plt.subplots(figsize=(7,5))
    types = ['dp_' + recoil_type + '/', 'dp_' + recoil_type + '/1bump/', 'dp_' + recoil_type + '/2bump/']
    hist_range = (0,100)
    for type in types:
        root_file = '/mnt/raid/users/duylai/warios_' + type + 'SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data = np.array(tree[histo_name], dtype=np.float64)
        plt.hist(data, bins=200, range=hist_range, histtype='step')
    ax.grid(alpha=0.5, zorder=0)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)
    plt.legend([r'All events', r'$\log(S2/S1)\geq0.5$, DP', r'$\log(S2/S1)<0.5$, DP'], loc=0, fontsize=font_size)
    plt.xlim(0,120)
    # plt.ylim(0,1e6)
    plt.title(r'Inverse RoI', fontsize=font_size, color='r')
    plt.xlabel(r'S1 [photoelectrons]', fontsize=font_size+2)
    plt.ylabel(r'Counts', fontsize=font_size+2)
    fig.tight_layout()
    if save_fig:
        plt.savefig('python/' + histo_name + '_inverse_RoI.png', dpi=300)
    if show_fig:
        plt.show()

elif histo_name[4:] == 'F90_vs_S1L':
    fig, ax = plt.subplots(figsize=(7,5))
    root_file = '/mnt/raid/users/duylai/warios_' + type + '/SuperWario.root'
    histo = uproot.open(root_file)[histo_name]
    data = histo.to_numpy()
    X, Y = np.meshgrid(data[2], data[1])
    Z = data[0]
    plt.pcolor(Y, X, Z, norm=colors.LogNorm(), cmap='jet')
    if histo_name[:3] == 'RoI':
        pre_histo = uproot.open(root_file)['pre_F90_vs_S1L']
        plt.title(r'RoI cut efficiency = ' + str(round(histo.member('fEntries')/pre_histo.member('fEntries')*100, 5))
                                        + r' %', fontsize=font_size-4)
    ax.grid(alpha=0.5, zorder=0)
    cb = plt.colorbar()
    cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
    plt.xlabel(r'$S1$ [photoelectrons]', fontsize=font_size+2)
    plt.ylabel(r'$F90$', fontsize=font_size+2)
    if histo_name[:3] == 'pre':
        plt.xlim(0,800)
        # RoI_cuts()
    elif histo_name[:3] == 'RoI':
        plt.xlim(0,125)
    plt.ylim(0,1)
    fig.tight_layout()
    if save_fig:
        plt.savefig('python/' + histo_name + '_' + type + '.png', dpi=300)
    if show_fig:
        plt.show()

elif histo_name[4:] == 'logS2maxLoS1L_vs_F90':
    fig, ax = plt.subplots(figsize=(7,5))
    root_file = '/mnt/raid/users/duylai/warios_' + type + '/SuperWario.root'
    histo = uproot.open(root_file)[histo_name]
    data = histo.to_numpy()
    X, Y = np.meshgrid(data[2], data[1])
    Z = data[0]
    plt.pcolor(Y, X, Z, norm=colors.LogNorm(), cmap='jet')
    if histo_name[:3] == 'RoI':
        pre_histo = uproot.open(root_file)['pre_logS2maxLoS1L_vs_F90']
        plt.title(r'RoI cut efficiency = ' + str(round(histo.member('fEntries')/pre_histo.member('fEntries')*100, 5))
                                        + r' %', fontsize=font_size-4)
    ax.grid(alpha=0.5, zorder=0)
    cb = plt.colorbar()
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    plt.xlabel(r'$F90$', fontsize=font_size+2)
    plt.ylabel(r'$\log_{10}\left(\frac{S2}{S1}\right)$', fontsize=font_size+2)
    plt.xlim(0,1)
    plt.ylim(-3,5)
    fig.tight_layout()
    if save_fig:
        plt.savefig('python/' + histo_name + '_' + type + '.png', dpi=300)
    if show_fig:
        plt.show()