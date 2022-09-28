#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import (MultipleLocator, LogLocator, ScalarFormatter)
import uproot

save_fig = 1

# recoil_type = 'ER'
# recoil_type = 'NR'
recoil_type = ''
# data_type = 'dp'
data_type = 'neutron'
# cut = '++'
cut = ''

# names = ['S2maxDT', 'S2maxLY', 'S2maxFracS2L', 'S2maxQMatch', 'S2maxLoS1Lcorr_vs_S1F90', 'S1TTR', 'S1TTT',
            # 'S1maxFracPMT', 'S1LY', 'logS2maxLoS1L_vs_QMatch', 'S1Nx_vs_S2maxFracS2L', 'S2maxLSY_vs_S2maxLSX']
# names = ['nS2_vs_S2maxFracS2L', 'S1Nx_vs_S2maxFracS2L', 'nS2_vs_S2maxFracS2L_ratio',
            # 'S1Nx_vs_S2maxFracS2L_ratio', 'S2maxLoS1Lcorr_vs_S1F90']
names = ['S2maxLoS1Lcorr_vs_S1F90', 'S2maxDT', 'S2maxDT_vs_Radius', 'S2maxLoS1Lcorr_vs_S2maxDT']

if save_fig == 0:
    show_fig = 1
else:
    show_fig = 0

# Plotting fit_paramters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size + 2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size - 4
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams["legend.markerscale"] = 3.5
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)

def cut_on_ratio(bump, var_x, var_y, var_to_cut_x, var_to_cut_y, thresh):
    if bump == 'both':
        data_cut_upper_x, data_cut_upper_y = cut_on_ratio('upper', var_x, var_y, var_to_cut_x, var_to_cut_y, thresh)
        data_cut_lower_x, data_cut_lower_y = cut_on_ratio('lower', var_x, var_y, var_to_cut_x, var_to_cut_y, thresh)
        return np.append(data_cut_upper_x, data_cut_lower_x), np.append(data_cut_upper_y, data_cut_lower_y)
    #################### upper bump ####################
    root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
    if cut == '++':
        root_file += '_++'
    root_file += '/1bump/SuperWario.root'
    tree = uproot.open(root_file + ':HebingTree')
    if bump == 'upper':
        data_x = tree[var_x].array(library='np')
        data_y = tree[var_y].array(library='np')
    data_to_cut_upper_x = tree[var_to_cut_x].array(library='np')
    data_to_cut_upper_y = tree[var_to_cut_y].array(library='np')
    #################### lower bump ####################
    root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
    if cut == '++':
        root_file += '_++'
    root_file += '/2bump/SuperWario.root'
    tree = uproot.open(root_file + ':HebingTree')
    if bump == 'lower':
        data_x = tree[var_x].array(library='np')
        data_y = tree[var_y].array(library='np')
    data_to_cut_lower_x = tree[var_to_cut_x].array(library='np')
    data_to_cut_lower_y = tree[var_to_cut_y].array(library='np')
    #################### selecting by ratio ####################
    val_upper, edges_x , edges_y = np.histogram2d(data_to_cut_upper_x, data_to_cut_upper_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
    val_lower, _ , _ = np.histogram2d(data_to_cut_lower_x, data_to_cut_lower_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
    ratio = np.nan_to_num(np.divide(val_lower, val_upper))
    data_cut_x = np.array([])
    data_cut_y = np.array([])
    for i in np.arange(len(data_x)):
        if bump == 'upper':
            ind_x = np.digitize(data_to_cut_upper_x[i], edges_x, right=True)
            ind_y = np.digitize(data_to_cut_upper_y[i], edges_y, right=True)
        elif bump == 'lower':
            ind_x = np.digitize(data_to_cut_lower_x[i], edges_x, right=True)
            ind_y = np.digitize(data_to_cut_lower_y[i], edges_y, right=True)
        if ratio[ind_x-1, ind_y-1] >= thresh:
            data_cut_x = np.append(data_cut_x, data_x[i])
            data_cut_y = np.append(data_cut_y, data_y[i])
    return data_cut_x, data_cut_y

# x, y = cut_on_ratio('lower', 'S1F90', 'S2maxLoS1Lcorr', 'S2maxFracS2L', 'nS2', 2.2)
# print(x.shape)
# print(y.shape)

for histo_name in names:
    if histo_name == 'S2maxDT':
        fig, ax = plt.subplots(figsize=(7,5))
        # bumps = ['', '/1bump', '/2bump']
        data_folder = ['neutron']
        nentries = np.array([], dtype=np.int)
        for data_type in data_folder:
        # for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type # + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=100, range=(0,1.4), histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            # if recoil_type == 'ER':
                # plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                # plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            # elif recoil_type == 'NR':
                # plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                # plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            plt.legend([r'Neutron runs, N = {:.2e}'.format(nentries[0])])
            plt.savefig('figures_neutron/' + histo_name + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S2LY' or histo_name == 'S2maxLY':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=200, range=(0,2000), histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        ax.set_yscale('log')
        # formatter = ScalarFormatter(useMathText=True)
        # formatter.set_scientific(True) 
        # formatter.set_powerlimits((-1,1))
        # ax.yaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S2maxFracS2L':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=200, range=(0.5,1), log=1, histtype='step')
        if recoil_type == 'ER':
            plt.ylim(1e3,1e8)
        ax.grid(alpha=0.5, zorder=0)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S2maxQMatch':
        bumps = ['', '/1bump', '/2bump']
        S1Nx = [0,1,2,3,4,5]
        nS2 = [1,2,3,4,5]
        from matplotlib import cm
        from scipy.optimize import curve_fit
        def gauss(x, a, c):
            return a*np.exp(-0.5*((x-1.0)/c)**2)
        for bump in bumps:
            fig, ax = plt.subplots(figsize=(7,5))
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            cut_var = tree['S1Nx'].array(library='np')
            c = np.arange(1, len(S1Nx) + 1)
            cnorm = colors.Normalize(vmin=c.min(), vmax=c.max())
            cmap = cm.ScalarMappable(norm=cnorm, cmap=cm.Blues)
            cmap.set_array([])
            for i in S1Nx:
                data_cut = data[np.where(cut_var == i)]
                if bump == '/2bump' and cut == '':
                    bin_val, edges = np.histogram(data_cut, bins=100, range=(0,4))
                    bin_pos = (edges[:-1]+edges[1:])/2
                    fit_param, _ = curve_fit(gauss, bin_pos, bin_val)
                    fwmh = 2.355*fit_param[1]
                    plt.hist(data_cut, bins=100, range=(0,4), histtype='step',
                        label='S1Nx = ' + str(i)  + ',    FWHM = ' + '{:.3f}'.format(np.abs(2.3548200450309493*fit_param[1])),
                        color=cmap.to_rgba(len(S1Nx) - i))
                    plt.plot(bin_pos, gauss(bin_pos, fit_param[0], fit_param[1]), linestyle=':', c='g', alpha=0.5)
                else:
                    plt.hist(data_cut, bins=100, range=(0,4), histtype='step', label='S1Nx = ' + str(i),
                        color=cmap.to_rgba(len(S1Nx) - i))
            cut_var = tree['nS2'].array(library='np')
            c = np.arange(1, len(nS2) + 1)
            cnorm = colors.Normalize(vmin=c.min(), vmax=c.max())
            cmap = cm.ScalarMappable(norm=cnorm, cmap=cm.Reds)
            cmap.set_array([])
            for i in nS2:
                data_cut = data[np.where(cut_var == i)]
                if bump == '/2bump' and cut == '':
                    bin_val, edges = np.histogram(data_cut, bins=100, range=(0,4))
                    bin_pos = (edges[:-1]+edges[1:])/2
                    fit_param, _ = curve_fit(gauss, bin_pos, bin_val)
                    plt.hist(data_cut, bins=100, range=(0,4), histtype='step',
                        label='nS2 = ' + str(i) + ',    FWHM = ' + '{:.3f}'.format(np.abs(2.3548200450309493*fit_param[1])),
                        color=cmap.to_rgba(len(nS2) - i))
                    plt.plot(bin_pos, gauss(bin_pos, fit_param[0], fit_param[1]), linestyle=':', c='g', alpha=0.5)
                else: 
                    plt.hist(data_cut, bins=100, range=(0,4), histtype='step', label='nS2 = ' + str(i),
                        color=cmap.to_rgba(len(nS2) - i))
            ax.grid(alpha=0.5, zorder=0)
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1))
            ax.yaxis.set_major_formatter(formatter)
            plt.legend()
            plt.title(cut + r', ' + data_type, fontsize=font_size,)
            plt.xlabel(histo_name, fontsize=font_size+2, color='r')
            plt.ylabel(r'Counts', fontsize=font_size+2)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '.png', dpi=300)
                elif recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S2maxLoS1Lcorr_vs_S1F90':
        bumps = ['']
        # bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type # + '_' + recoil_type
            if cut == '++':
                pre_root_file = root_file
                root_file += '_++'
                pre_root_file += bump + '/SuperWario.root'
                pre_tree = uproot.open(pre_root_file + ':HebingTree')
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S1F90'].array(library='np')
            data_y = tree['S2maxLoS1Lcorr'].array(library='np')
            thresh = 3,
            # if bump == '':
                # data_x, data_y = cut_on_ratio('both', 'S1F90', 'S2maxLoS1Lcorr', 'S2maxFracS2L', 'nS2', thresh)
            # elif bump == '/1bump':
                # data_x, data_y = cut_on_ratio('upper', 'S1F90', 'S2maxLoS1Lcorr', 'S2maxFracS2L', 'nS2', thresh)
            # elif bump == '/2bump':
                # data_x, data_y = cut_on_ratio('lower', 'S1F90', 'S2maxLoS1Lcorr', 'S2maxFracS2L', 'nS2', thresh)
            if cut == '++':
                print(r'Cut efficiency ', bump[1:], r': ', 100*len(data_y)/len(pre_tree['S2maxLoS1Lcorr'].array(library='np')))
            # drifttime = tree['S2maxDT'].array(library='np')
            # driftcorr = np.exp(-(drifttime)/0.5)
            # data_y = np.divide(data_y,driftcorr)
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[200,200], range=[[0,1],[-3,5]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            plt.title(r'{}, {}, N = {:.2e}'.format(cut, data_type, data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S1F90', fontsize=font_size+2)
            plt.ylabel(r'$\log_{10}\left(\frac{S2}{S1}\right)$', fontsize=font_size+2)
            fig.tight_layout()
            if save_fig:
                # if recoil_type == 'ER':
                    # plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                # elif recoil_type == 'NR':
                    # plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                plt.savefig('figures_' + data_type + '/' + histo_name + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1maxFracPMT':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=100, range=(0,0.4), histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S1TTR':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=100, range=(0,1), histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S1TTT':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=100, histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S1LY':
        fig, ax = plt.subplots(figsize=(7,5))
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data = tree[histo_name].array(library='np')
            nentries = np.append(nentries, data.shape[0])
            plt.hist(data, bins=100, range=(0,100), histtype='step')
        ax.grid(alpha=0.5, zorder=0)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(histo_name, fontsize=font_size+2, color='r')
        plt.ylabel(r'Counts', fontsize=font_size+2)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.5$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.5$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            elif recoil_type == 'NR':
                plt.legend([r'All events, N = {:.2e}'.format(nentries[0]), r'$\log(S2/S1)\geq0.2$, N = {:.2e}'.format(nentries[1]), r'$\log(S2/S1)<0.2$, N = {:.2e}'.format(nentries[2])], loc=0, fontsize=font_size)
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()
    
    elif histo_name == 'S2maxLSY_vs_S2maxLSX':
        bumps = ['', '/1bump', '/2bump']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree[histo_name].array(library='np')
            data_y = tree['S2maxLSY'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[200,200], range=[[-500,500],[-500,500]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxLSX', fontsize=font_size+2)
            plt.ylabel(r'S2maxLSY', fontsize=font_size+2)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'logS2maxLoS1L_vs_QMatch':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxQMatch'].array(library='np')
            data_y = tree['S2maxLoS1Lcorr'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[200,200], range=[[0.,4.],[-3.,5.],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.yaxis.set_major_locator(MultipleLocator(1.0))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxQMatch', fontsize=font_size+2)
            plt.ylabel(r'$\log_{10}\left(\frac{S2}{S1}\right)$', fontsize=font_size+2)
            plt.xlim(0,4)
            plt.ylim(-3,5)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1Nx_vs_S2maxQMatch':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxQMatch'].array(library='np')
            data_y = tree['S1Nx'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[100,11], range=[[0.,4.],[-0.5,10.5],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.yaxis.set_major_locator(MultipleLocator(1))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxQMatch', fontsize=font_size+2)
            plt.ylabel(r'S1Nx', fontsize=font_size+2)
            plt.xlim(0.,4.)
            plt.ylim(-1,11)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1Nx_vs_S2maxFracS2L':
        bumps = ['', '/1bump', '/2bump']
        # ratio = compute_ratio('S2maxFracS2L', 'S1Nx')
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxFracS2L'].array(library='np')
            data_y = tree['S1Nx'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            # val, edges_x, edges_y = np.histogram2d(data_x, data_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
            # for i in np.arange(100):
                # for j in np.arange(11):
                    # if ratio[i,j] < 1.5:
                        # val[i,j] = 0
            # X, Y = np.meshgrid(edges_y[:], edges_x[:])
            # plt.pcolor(Y, X, val, norm=colors.LogNorm(), cmap='jet')
            plt.hist2d(data_x, data_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxFracS2L', fontsize=font_size+2)
            plt.ylabel(r'S1Nx', fontsize=font_size+2)
            plt.xlim(0.,1.)
            plt.ylim(-1,11)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1Nx_vs_S2maxFracS2L_ratio':
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/1bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_upper_x = tree['S2maxFracS2L'].array(library='np')
        data_upper_y = tree['S1Nx'].array(library='np')
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/2bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_lower_x = tree['S2maxFracS2L'].array(library='np')
        data_lower_y = tree['S1Nx'].array(library='np')
        fig, ax = plt.subplots(figsize=(7,5))
        val_upper, edges_x, edges_y = np.histogram2d(data_upper_x, data_upper_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
        val_lower, _ , _ = np.histogram2d(data_lower_x, data_lower_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
        X, Y = np.meshgrid(edges_y[:], edges_x[:])
        plt.pcolor(Y, X, np.divide(val_lower,val_upper), cmap='Blues', vmax=3.5)
        ax.grid(alpha=0.5, zorder=0)
        cb = plt.colorbar()
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(r'S2maxFracS2L', fontsize=font_size+2)
        plt.ylabel(r'S1Nx', fontsize=font_size+2)
        plt.xlim(0.,1.)
        plt.ylim(-1,11)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            if recoil_type == 'NR':
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'S2maxFracS2L_vs_S2maxQMatch':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxQMatch'].array(library='np')
            data_y = tree['S2maxFracS2L'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[100,100], range=[[0.,4.],[0.,1.],], norm=colors.LogNorm(), cmap='jet')
            ax.yaxis.set_major_locator(MultipleLocator(0.1))
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxQMatch', fontsize=font_size+2)
            plt.ylabel(r'S2maxFracS2L', fontsize=font_size+2)
            plt.xlim(0.,4.)
            plt.ylim(0.,1.)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1Nx_vs_nS2':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['nS2'].array(library='np')
            data_y = tree['S1Nx'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[11,11], range=[[-0.5,10.5],[-0.5,10.5],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'nS2', fontsize=font_size+2)
            plt.ylabel(r'S1Nx', fontsize=font_size+2)
            plt.xlim(-1,11)
            plt.ylim(-1,11)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'S1Nx_vs_nS2_ratio':
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/1bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_upper_x = tree['nS2'].array(library='np')
        data_upper_y = tree['S1Nx'].array(library='np')
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/2bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_lower_x = tree['nS2'].array(library='np')
        data_lower_y = tree['S1Nx'].array(library='np')
        fig, ax = plt.subplots(figsize=(7,5))
        val_upper, edges_x, edges_y = np.histogram2d(data_upper_x, data_upper_y, bins=[11,11], range=[[-0.5,10.5],[-0.5,10.5],])
        val_lower, _ , _ = np.histogram2d(data_lower_x, data_lower_y, bins=[11,11], range=[[-0.5,10.5],[-0.5,10.5],])
        X, Y = np.meshgrid(edges_y[:], edges_x[:])
        plt.pcolor(Y, X, np.divide(val_lower,val_upper), cmap='Blues', vmax=3.5)
        ax.grid(alpha=0.5, zorder=0)
        cb = plt.colorbar()
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(r'nS2', fontsize=font_size+2)
        plt.ylabel(r'S1Nx', fontsize=font_size+2)
        plt.xlim(-1,11)
        plt.ylim(-1,11)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            if recoil_type == 'NR':
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'nS2_vs_S2maxFracS2L':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxFracS2L'].array(library='np')
            data_y = tree['nS2'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxFracS2L', fontsize=font_size+2)
            plt.ylabel(r'nS2', fontsize=font_size+2)
            plt.xlim(0.,1.)
            plt.ylim(-1,11)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()

    elif histo_name == 'nS2_vs_S2maxFracS2L_ratio':
        # bumps = ['', '/1bump', '/2bump']
        # for bump in bumps:
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/1bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_upper_x = tree['S2maxFracS2L'].array(library='np')
        data_upper_y = tree['nS2'].array(library='np')
        root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
        if cut == '++':
            root_file += '_++'
        root_file += '/2bump/SuperWario.root'
        tree = uproot.open(root_file + ':HebingTree')
        data_lower_x = tree['S2maxFracS2L'].array(library='np')
        data_lower_y = tree['nS2'].array(library='np')
        fig, ax = plt.subplots(figsize=(7,5))
        val_upper, edges_x, edges_y = np.histogram2d(data_upper_x, data_upper_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
        val_lower, _ , _ = np.histogram2d(data_lower_x, data_lower_y, bins=[100,11], range=[[0.,1.],[-0.5,10.5],])
        X, Y = np.meshgrid(edges_y[:], edges_x[:])
        plt.pcolor(Y, X, np.divide(val_lower,val_upper), cmap='Blues', vmax=3.)
        ax.grid(alpha=0.5, zorder=0)
        cb = plt.colorbar()
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        plt.title(cut + r', ' + data_type, fontsize=font_size,)
        plt.xlabel(r'S2maxFracS2L', fontsize=font_size+2)
        plt.ylabel(r'nS2', fontsize=font_size+2)
        plt.xlim(0.,1.)
        plt.ylim(-1,11)
        fig.tight_layout()
        if save_fig:
            if recoil_type == 'ER':
                plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + cut + '.png', dpi=300)
            if recoil_type == 'NR':
                plt.savefig('pileup/RoI/' + histo_name + '_' + cut + '.png', dpi=300)
        if show_fig:
            plt.show()

    elif histo_name == 'nS2_vs_S2maxQMatch':
        bumps = ['', '/1bump', '/2bump']
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxQMatch'].array(library='np')
            data_y = tree['nS2'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[100,11], range=[[0.,4.],[-0.5,10.5],], norm=colors.LogNorm(), cmap='jet')
            ax.grid(alpha=0.5, zorder=0)
            cb = plt.colorbar()
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxQMatch', fontsize=font_size+2)
            plt.ylabel(r'nS2', fontsize=font_size+2)
            plt.xlim(0.,4.)
            plt.ylim(-1,11)
            fig.tight_layout()
            if save_fig:
                if recoil_type == 'ER':
                    plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                if recoil_type == 'NR':
                    plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
            if show_fig:
                plt.show()
    
    elif histo_name == 'S2maxDT_vs_Radius':
        bumps = ['']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type # + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = np.sqrt(tree['S2maxLSX'].array(library='np')**2+tree['S2maxLSY'].array(library='np')**2)
            data_y = tree['S2maxDT'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[200,200], range=[[0,500],[0,1.4]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'Radius', fontsize=font_size+2)
            plt.ylabel(r'S2maxDT', fontsize=font_size+2)
            fig.tight_layout()
            if save_fig:
                # if recoil_type == 'ER':
                    # plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                # elif recoil_type == 'NR':
                    # plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                plt.savefig('figures_' + data_type + '/' + histo_name + '.png', dpi=300)
            if show_fig:
                plt.show()
    
    elif histo_name == 'S2maxLoS1Lcorr_vs_S2maxDT':
        bumps = ['']
        nentries = np.array([], dtype=np.int)
        for bump in bumps:
            root_file = '/mnt/raid/users/duylai/warios_' + data_type # + '_' + recoil_type
            if cut == '++':
                root_file += '_++'
            root_file += bump + '/SuperWario.root'
            tree = uproot.open(root_file + ':HebingTree')
            data_x = tree['S2maxDT'].array(library='np')
            data_y = tree['S2maxLoS1Lcorr'].array(library='np')
            fig, ax = plt.subplots(figsize=(7,5))
            plt.hist2d(data_x, data_y, bins=[200,200], range=[[0,1.2],[-3,5]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            plt.title(cut + r', ' + data_type + r', N = {:.2e}'.format(data_x.shape[0]), fontsize=font_size,)
            plt.xlabel(r'S2maxDT', fontsize=font_size+2)
            plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+2)
            fig.tight_layout()
            if save_fig:
                # if recoil_type == 'ER':
                    # plt.savefig('pileup/inverse_RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                # elif recoil_type == 'NR':
                    # plt.savefig('pileup/RoI/' + histo_name + '_' + bump[1:] + '_' + cut + '.png', dpi=300)
                plt.savefig('figures_' + data_type + '/' + histo_name + '.png', dpi=300)
            if show_fig:
                plt.show()