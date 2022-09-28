######### Taking any .root file with a HebingTree as input data, then draw the plots #########
# written by Duy Lai, 09.2022

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors, patches
from matplotlib.ticker import (MultipleLocator, LogLocator, ScalarFormatter)
import uproot, os, h5py

plot_thesis_1col = 0
plot_thesis_2col = 1

# save_fig is a boolean, if False, show_plot will be True and vice versa
save_fig = 1
# save as .png or .svg (the latter for vector graphics, not recommended for the 2D plots with many colors)
file_type = '.png'

nbins = 200

names = []
# names = ['S2maxLY_vs_RunDate']
# names = ['S2maxLoS1Lcorr_vs_RunDate']
# names = ['S1TTR_vs_S2maxDT']
# names = ['S2maxLoS1Lcorr_vs_S1F90_runs']
# names = ['S1LY', 'S1TTR']
# names = ['eff_detec_rec']
# names = ['e_depo']
names = ['eff_match']
names = ['S2maxDT_vs_Radius']
# names = ['S2maxLoS1Lcorr_vs_S1F90_S1bins']
# names = ['wario_efficiency_S1bins']
# names = ['S1TTR_vs_S1F90', 'S1TTR_vs_S2maxDT', 'S2maxQMatch_vs_S2maxDT']
# names = ['S1TTR', 'S1maxFracPMT', 'S1LY']
# names = ['nS2', 'nS1perS2max', 'S2maxQMatch']
# names = ['S2maxLSY_vs_S2maxLSX']
# names = ['S2maxLSY_vs_S2maxLSX', 'S2maxDT', 'S2maxDT_vs_Radius', 'S2maxLoS1Lcorr_vs_S2maxDT']
# names = ['S1LY', 'S1TTR', 'S1F90_vs_S1LY', 'S2maxLoS1Lcorr_vs_S1F90', 'S2maxLSY_vs_S2maxLSX', 'S2maxDT',
    #   'S2maxDT_vs_Radius', 'S2maxLoS1Lcorr_vs_S2maxDT']

# List of data types/sources
data_folder = ['dp_earlies', 'dp_latest']
# data_folder = ['neutron']
# data_folder = ['DP2']

input_path = '/mnt/raid/users/duylai/warios_'
# input_path = '/mnt/raid/users/duylai/preselection_'
# input_path = '/mnt/raid/users/duylai/matched_'

if plot_thesis_1col or plot_thesis_2col:
    output_path = 'thesis_plots/'
else:
    output_path = 'warios_plots/'

# For each type of data, create a DataFrame from the SuperWario then group them into a collection/dictionary
def import_data(Data_Collection, folder, data_type):
    df = pd.DataFrame()
    tree = uproot.open(folder + 'SuperWario.root:HebingTree')
    for var in tree.keys():
        df.insert(df.shape[1], var, tree[var].array(library='np'))
        # print(r'{} -- {}'.format(var, tree[var].array(library='np').mean()))
    Data_Collection[data_type] = df

Data_Collection = {}
for data_type in data_folder:
    import_data(Data_Collection, input_path + data_type + '/', data_type)

# Create a csv file for the runs' time and ID, then import it as a DataFrame
csv_name = 'All'
All_runs_csv_path = csv_name + '_runs.csv'
if os.path.exists(All_runs_csv_path):
    All_runs = pd.read_csv(All_runs_csv_path)
else:
    fh5 = h5py.File(csv_name + '_runs.h5', 'r')
    fh5_runIDs = np.array(fh5['runID'][:], dtype=int)
    fh5_DAQTimeLapse = np.array(fh5['DAQTimeLapse'][:], dtype=np.float64)
    fh5_Drift = np.array(fh5['DriftField_V_per_cm'][:], dtype=np.float64)
    fh5_dates = np.array([], dtype = np.datetime64)
    for t in fh5['DAQStartTime'][:]:
        t = t.decode('utf-8')
        dt = np.array([t[6:10] + '-' + t[0:2] + '-' + t[3:5] + 'T' + t[12:20]], dtype=np.datetime64)
        fh5_dates = np.append(fh5_dates, dt)
    pd.DataFrame({'RunID': fh5_runIDs, 'Date': fh5_dates, 'DAQTimeLapse': fh5_DAQTimeLapse,
            'DriftField_V_per_cm': fh5_Drift}).to_csv(All_runs_csv_path, index=False)

fig_size = (7,5)
if plot_thesis_1col:
    fig_size = (8,5)
    mag = -1
elif plot_thesis_2col:
    mag = 6
else:
    mag = 2

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
plt.rc('xtick', labelsize=font_size + mag)
plt.rc('ytick', labelsize=font_size + mag)

def identify_dates(runID, unit='s'):
    return np.datetime64(All_runs.query('RunID == {}'.format(runID))['Date'].to_numpy()[0], unit)

# Formatting floats for StatsBox
def format_number(num):
    x = int(np.log10(abs(num)))
    if num < 1:
        x -= 1
    if x <= -3:
        return r'{:>.2e}'.format(num)
    else:
        return (r'{:>.' + str(3-x) + r'f}').format(num)

# Produce the ROOT-like statistical box
def StatsBox(dd = True):
    legend_size = font_size
    if plot_thesis_2col:
        legend_size = font_size + mag -2
    if dd == True:
        if not plot_thesis_2col or plot_name == 'S2maxLoS1Lcorr_vs_S1F90':
            handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 5
            labels = []
            labels.append(r'Entries           {:,}'.format(data_x.shape[0]))
            labels.append(r'Mean x           ' + format_number(data_x.astype('float64').mean()))
            labels.append(r'Mean y           ' + format_number(data_y.astype('float64').mean()))
            labels.append(r'Std Dev x       ' + format_number(data_x.astype('float64').std()))
            labels.append(r'Std Dev y       ' + format_number(data_y.astype('float64').std()))
            ax.legend(handles, labels, loc='best', fontsize= legend_size, 
                fancybox=True, framealpha=0.5, 
                handlelength=0, handletextpad=0)
    else:
        handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 5
        labels = []
        labels.append(r'Entries         {:,}'.format(data.shape[0]))
        labels.append(r'Mean           ' + format_number(data.astype('float64').mean()))
        labels.append(r'Std Dev       ' + format_number(data.astype('float64').std()))
        loca = 'best'
        if plot_thesis_2col:
            if plot_name == 'S1LY':
                loca = 'upper left'
            elif plot_name == 'S1TTR':
                loca = 'lower right'
        ax.legend(handles, labels, loc=loca, fontsize= legend_size, 
            fancybox=True, framealpha=0.5, 
            handlelength=0, handletextpad=0)

def Title():
    if not plot_thesis_1col and not plot_thesis_2col:
        if data_type == 'neutron':
            plt.title(r'DP2 Neutron runs', fontsize=font_size, color='r')
        elif data_type == 'DP1':
            plt.title(r'All DP1 runs', fontsize=font_size, color='r')
        elif data_type == 'DP2':
            plt.title(r'All DP2 runs', fontsize=font_size, color='r')

def f_indicate_date_ranges():
    locator = mdates.DayLocator(bymonthday=(1, 4, 7, 10, 13, 16, 19, 22, 25, 28))
    # locator = mdates.DayLocator(bymonthday=(1, 5, 9, 13, 17, 21, 25, 29))
    formatter = mdates.ConciseDateFormatter(locator, show_offset=False)
    alpha = 0.3
    y_limits = plt.ylim()
    plt.xlim(np.datetime64('2019-06-01'), np.datetime64('2019-08-28'))
    date_proxy = np.arange(np.datetime64('2019-05-12T14:35:49'), np.datetime64('2019-08-26T12:51:30'),
                np.timedelta64(30, 's'))
    # SP
    plt.fill_between(x=date_proxy, y1=y_limits[0], y2=y_limits[1],
                     where=(np.datetime64('2019-05-12T14:35:49') <= date_proxy)
                           & (date_proxy <= np.datetime64('2019-06-05T11:54:21')),
                     facecolor='green', edgecolor='none',  alpha=alpha, label='SP')
    plt.fill_between(x=date_proxy, y1=y_limits[0], y2=y_limits[1],
                     where=(np.datetime64('2019-08-26T01:40:57') <= date_proxy)
                           & (date_proxy <= np.datetime64('2019-08-26T12:51:30')),
                     facecolor='green', edgecolor='none', alpha=alpha)

    # SPD
    plt.fill_between(x=date_proxy, y1=y_limits[0], y2=y_limits[1],
                     where=(np.datetime64('2019-06-19T10:49:18') <= date_proxy)
                           & (date_proxy <= np.datetime64('2019-06-19T18:29:03')),
                     facecolor='blue', edgecolor='none', alpha=alpha, label='SPD')

    # HV failure
    plt.plot([np.datetime64('2019-07-07T15:55:13'), np.datetime64('2019-07-07T15:55:13')], y_limits, color='black',
             linewidth=1.0, alpha=alpha, label='HV failure')

    # DAQ Problems
    plt.fill_between(x=date_proxy, y1=y_limits[0], y2=y_limits[1],
                     where=(np.datetime64('2019-08-01T11:04:36') <= date_proxy) &
                           (date_proxy <= np.datetime64('2019-08-05T11:05:23')),
                     facecolor='red', edgecolor='none', alpha=alpha, label='DAQ problems')
    plt.xticks(rotation=50)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    return

if save_fig == 0:
    show_fig = 1
else:
    show_fig = 0

for plot_name in names:

    if plot_name in ['S2maxDT', 'S1LY', 'S2maxQMatch', 'S1TTR', 'S1maxFracPMT', 'nS2', 'nS1perS2max', 'S2maxQMatch']:
        # data_folder = ['neutron']
        # data_folder = ['neutron', 'dp_earlies', 'dp_latest']
        if plot_name == 'S2maxDT':
            xrange = (0,1.2)
        elif plot_name == 'S1LY':
            xrange = (0,100)
            # xrange = (0,800)
        elif plot_name == 'S2maxQMatch':
            xrange = (0.0,4.0)
        elif plot_name == 'S1TTR':
            xrange = (0.0,1.0)
        elif plot_name == 'S1maxFracPMT':
            xrange = (0.0,1.0)
        elif plot_name == 'nS2':
            xrange = (0,10)
        elif plot_name == 'nS1perS2max':
            xrange = (0,10)
        elif plot_name == 'S2maxQMatch':
            xrange = (0.0,4.0)
        fig, ax = plt.subplots(figsize=fig_size)
        for data_type in data_folder:
            if (data_type not in Data_Collection.keys()):
                import_data(Data_Collection, input_path + data_type + '/', data_type)
            data = Data_Collection[data_type][plot_name]
            if data_type == 'neutron':
                label = r'DP2 Neutron runs'
            elif data_type == 'DP1':
                label = r'All DP1 runs'
            elif data_type == 'DP2':
                label = r'All DP2 runs'
            elif data_type == 'dp_earlies':
                label = r'Early DP runs'
            elif data_type == 'dp_latest':
                label = r'Late DP runs'
            if plot_name in ['nS2', 'nS1perS2max']:
                bins = plt.hist(data, bins=10, align='left', range=xrange, density=False, histtype='step', label=label)
            else:
                bins = plt.hist(data, bins=100, range=xrange, density=False, histtype='step', label=label)
            print(data.mean())
            print(data.std())
            # bins = plt.hist(data, bins=100, range=xrange, density=True, histtype='step', label=label)[1]
        ax.grid(alpha=0.5, zorder=0)
        plt.xlim(xrange)
        yticks = plt.yticks()[0]
        if (plt.ylim()[1] != yticks[-1]):
            plt.ylim(yticks[0], 2*yticks[-2]-yticks[-3])
        # ax.yaxis.set_major_locator(MultipleLocator(0.4))
        # formatter = ScalarFormatter(useMathText=True)
        # formatter.set_scientific(True)
        # formatter.set_powerlimits((-1,1))
        # ax.yaxis.set_major_formatter(formatter)
        # plt.legend(fontsize=font_size-2)
        loca = 'best'
        if plot_name == 'S1LY':
            loca = 'upper left' 
            plt.legend(['DP1, mean 70.65, std 21.14', 'DP2, mean 68.15, std 20.81'], fontsize=font_size + mag -2, loc=loca)
        elif plot_name == 'S1TTR':
            loca = 'lower right' 
            plt.legend(['DP1, mean 0.4030, std 0.1070', 'DP2, mean 0.3731, std 0.1143'], fontsize=font_size + mag -2, loc=loca)
        elif plot_name == 'S2maxDT':
            plt.legend(['DP1, mean 0.3470, std 0.2691', 'DP2, mean 0.4398, std 0.2924'], fontsize=font_size + mag -2,)
        # StatsBox(False)
        xl = plot_name
        if plot_name == 'S1LY':
            for k in np.linspace(10,90,9):
                plt.plot([k,k], plt.ylim(), '--', c='r')
            ax.xaxis.set_major_locator(MultipleLocator(10.0))
            # ax.xaxis.set_major_locator(MultipleLocator(100.0))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1))
            ax.yaxis.set_major_formatter(formatter)
            xl += r' [p.e.]'
        elif plot_name == 'S1TTR':
            # plt.plot([0.15,0.15], plt.ylim(), '--', c='r')
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1))
            ax.yaxis.set_major_formatter(formatter)
        elif plot_name == 'S1maxFracPMT':
            plt.yscale('log')
            plt.grid(which='major', ls='-')
            plt.grid(which='minor', ls='--', alpha=0.7)
            plt.ylim(1e0,1e8)
            plt.plot([0.3,0.3], plt.ylim(), '--', c='r')
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
        elif plot_name == 'nS2':
            # plt.plot([0.5,0.5], plt.ylim(), '--', c='r')
            plt.plot([1.5,1.5], plt.ylim(), '--', c='r')
            ax.xaxis.set_major_locator(MultipleLocator(1))
        elif plot_name == 'nS1perS2max':
            # plt.plot([0.5,0.5], plt.ylim(), '--', c='r')
            plt.plot([1.5,1.5], plt.ylim(), '--', c='r')
            ax.xaxis.set_major_locator(MultipleLocator(1))
        elif plot_name == 'S2maxQMatch':
            plt.yscale('log')
            plt.grid(which='major', ls='-')
            plt.grid(which='minor', ls='--', alpha=0.7)
            plt.plot([0.75,0.75], plt.ylim(), '--', c='r')
            plt.plot([1.25,1.25], plt.ylim(), '--', c='r')
            plt.ylim(1e0,1e6)
            # ind = np.where((bins[1] >= 0.75) & (bins[1] <= 1.25))
            # plt.fill_between(bins[1][ind], bins[0][ind], color='r', alpha=0.5)
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1))
            ax.yaxis.set_major_formatter(formatter)
        elif plot_name == 'S2maxDT':
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1))
            ax.yaxis.set_major_formatter(formatter)
            xl += r' [ms]'
        plt.xlabel(xl, fontsize=font_size+mag)
        plt.ylabel(r'Counts', fontsize=font_size+mag)
        # plt.ylabel(r'Probability / ({:.3f})'.format((bins[1][1]-bins[1][0])), fontsize=font_size+mag)
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()

    elif plot_name == 'S2maxLoS1Lcorr_vs_S1F90':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S1F90']
            data_y = Data_Collection[data_type]['S2maxLoS1Lcorr']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,1],[-3,5]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            plt.ylim(-3,5)
            ax.yaxis.set_major_locator(MultipleLocator(1.0))
            Title()
            plt.xlabel(r'S1F90', fontsize=font_size+mag)
            plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()

    elif plot_name == 'S2maxLoS1Lcorr_vs_S1F90_runs':
        data_folder = ['neutron', 'DP']
        for data_type in data_folder:
            runs = next(os.walk('/mnt/raid/users/duylai/warios_' + data_type + '/'), (None, None, []))[2]
            summary = pd.read_csv('/mnt/raid/users/duylai/Satoshi_' + data_type + '_runs/satoshi-summary.txt', sep='\t')
            for run in runs:
                try:
                    run_nbr = int(run[8:12])
                except ValueError:
                    continue
                root_file = '/mnt/raid/users/duylai/warios_' + data_type + '/' + run
                tree = uproot.open(root_file + ':HebingTree')
                data_x = tree['S1F90'].array(library='np')
                data_y = tree['S2maxLoS1Lcorr'].array(library='np')
                fig, ax = plt.subplots(figsize=fig_size)
                plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,1],[-3,5]], norm=colors.LogNorm(), cmap='jet')
                cb = plt.colorbar()
                cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
                ax.grid(alpha=0.5, zorder=0)
                plt.title(r'{}, RunID {} - {}'.format(run[5:-5], summary.query('SRunID  == ' + str(run_nbr))['StartRunID'].to_numpy()[0],
                            summary.query('SRunID  == ' + str(run_nbr))['EndRunID'].to_numpy()[0], fontsize=font_size))
                plt.xlabel(r'S1F90', fontsize=font_size+mag)
                plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+mag)
                StatsBox()
                fig.tight_layout()
                if save_fig:
                    plt.savefig(output_path + data_type + '_runs/' + str(run_nbr) + '.png', dpi=300)
                    plt.close()
                if show_fig:
                    plt.show()

    elif plot_name == 'S2maxLoS1Lcorr_vs_S1F90_S1bins':
        mag += 8
        S1_bins = np.linspace(0,100,11)
        for data_type in data_folder:
            for i in np.arange(S1_bins.shape[0]-1):
                Bin_Data = Data_Collection[data_type].query(r'{} <= S1LY < {}'.format(S1_bins[i], S1_bins[i+1]))
                data_x = Bin_Data['S1F90']
                data_y = Bin_Data['S2maxLoS1Lcorr']
                fig, ax = plt.subplots(figsize=fig_size)
                plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,1],[-3,5]], norm=colors.LogNorm(), cmap='jet')
                cb = plt.colorbar()
                cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
                ax.grid(alpha=0.5, zorder=0)
                plt.ylim(-3,5)
                # ax.yaxis.set_major_locator(MultipleLocator(1.0))
                plt.title(r'{}$\leq$S1LY<{} p.e.'.format(int(S1_bins[i]), int(S1_bins[i+1])), fontsize=font_size+mag-2)
                plt.xlabel(r'S1F90', fontsize=font_size+mag)
                plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+mag)
                StatsBox()
                fig.tight_layout()
                if save_fig:
                    plt.savefig(r'{}{}_S1bins/{}_{}.png'.format(output_path, data_type, int(S1_bins[i]), int(S1_bins[i+1])), dpi=300)
                    plt.close()
                if show_fig:
                    plt.show()

    elif plot_name == 'wario_efficiency_S1bins':
        S1_bins = np.arange(0,100,10)
        Satoshi = {}
        import_data(Satoshi, '/mnt/raid/users/duylai/Satoshi_neutron_runs/', 'neutron', 2)
        eff = np.array([])
        for i in np.arange(S1_bins.shape[0]-1):
            Wario_Data = Data_Collection['neutron'].query(r'{} <= S1LY < {}'.format(S1_bins[i], S1_bins[i+1]))
            Satoshi_Data = Satoshi['neutron'].query(r'{} <= S1LY < {}'.format(S1_bins[i], S1_bins[i+1]))
            # eff = np.append(eff, len(Wario_Data))
            eff = np.append(eff, len(Wario_Data)/len(Satoshi_Data))
        print(eff)
        fig, ax = plt.subplots(figsize=fig_size)
        plt.scatter(0.5*(S1_bins[:-1]+S1_bins[1:]), eff)
        plt.grid()
        plt.xlim(0,100)
        plt.ylim(0,0.2)
        ax.xaxis.set_major_locator(MultipleLocator(10.0))
        plt.xlabel('S1L [p.e]')
        plt.ylabel('Overall cut efficiency')
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()

    elif plot_name == 'S1F90_vs_S1LY':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S1LY']
            data_y = Data_Collection[data_type]['S1F90']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,100],[0,1]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            # ax.yaxis.set_major_locator(MultipleLocator(1.0))
            Title()
            plt.xlabel(r'S1L [p.e.]', fontsize=font_size+mag)
            plt.ylabel(r'S1F90', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S1TTR_vs_S1F90':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S1F90']
            data_y = Data_Collection[data_type]['S1TTR']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0.,1.],[0.,1.]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S1F90', fontsize=font_size+mag)
            plt.ylabel(r'S1TTR', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S1TTR_vs_S2maxDT':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S2maxDT']
            data_y = Data_Collection[data_type]['S1TTR']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0.,1.2],[0.,1.]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S2maxDT [ms]', fontsize=font_size+mag)
            plt.ylabel(r'S1TTR', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S2maxQMatch_vs_S2maxDT':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S2maxDT']
            data_y = Data_Collection[data_type]['S2maxQMatch']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0.,1.2],[0.,4.]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S2maxDT [ms]', fontsize=font_size+mag)
            plt.ylabel(r'S2maxQMatch', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S2maxDT_vs_Radius':
        for data_type in data_folder:
            data_x = np.sqrt(Data_Collection[data_type]['S2maxLSX']**2+Data_Collection[data_type]['S2maxLSY']**2)
            data_y = Data_Collection[data_type]['S2maxDT']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,500],[0,1.2]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S2maxLSRadius [mm]', fontsize=font_size+mag)
            plt.ylabel(r'S2maxDT [ms]', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S2maxLoS1Lcorr_vs_S2maxDT':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S2maxDT']
            data_y = Data_Collection[data_type]['S2maxLoS1Lcorr']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[0,1.2],[-3,5]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.yaxis.set_major_locator(MultipleLocator(1.0))
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S2maxDT [ms]', fontsize=font_size+mag)
            plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()
    
    elif plot_name == 'S2maxLSY_vs_S2maxLSX':
        for data_type in data_folder:
            data_x = Data_Collection[data_type]['S2maxLSX']
            data_y = Data_Collection[data_type]['S2maxLSY']
            fig, ax = plt.subplots(figsize=fig_size)
            plt.hist2d(data_x, data_y, bins=[nbins,nbins], range=[[-500,500],[-500,500]], norm=colors.LogNorm(), cmap='jet')
            cb = plt.colorbar()
            cb.ax.yaxis.set_major_locator(LogLocator(numticks=15))
            ax.xaxis.set_major_locator(MultipleLocator(250))
            ax.yaxis.set_major_locator(MultipleLocator(250))
            ax.grid(alpha=0.5, zorder=0)
            Title()
            plt.xlabel(r'S2maxLSX [mm]', fontsize=font_size+mag)
            plt.ylabel(r'S2maxLSY [mm]', fontsize=font_size+mag)
            StatsBox()
            fig.tight_layout()
            if save_fig:
                plt.savefig(output_path + plot_name + '_' + data_type + '.png', dpi=300)
                plt.close()
            if show_fig:
                plt.show()

    elif plot_name == 'S2maxLoS1Lcorr_vs_RunDate':
        data_folder = ['DP1', 'DP2', 'neutron']
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(alpha=0.5, zorder=0)
        for data_type in data_folder:
            data_x = np.array([], dtype=np.datetime64)
            mean_y = np.array([], dtype=np.float64)
            std_y = np.array([], dtype=np.float64)
            runs = next(os.walk('/mnt/raid/users/duylai/warios_' + data_type + '/'), (None, None, []))[2]
            summary = pd.read_csv('/mnt/raid/users/duylai/warios_' + data_type + '/satoshi-summary.txt', sep='\t')
            for run in runs:
                try:
                    run_nbr = int(run[8:12])
                except ValueError:
                    continue
                start_date = identify_dates(summary.query('SRunID  == ' + str(run_nbr))['StartRunID'].to_numpy()[0])
                end_date = identify_dates(summary.query('SRunID  == ' + str(run_nbr))['EndRunID'].to_numpy()[0])
                mid_date = np.datetime64(round(0.5 * (start_date.astype('int') + end_date.astype('int'))), 's')
                data_x = np.append(data_x, mid_date)
                root_file = '/mnt/raid/users/duylai/warios_' + data_type + '/' + run
                tree = uproot.open(root_file + ':HebingTree')
                data_y = tree['S2maxLoS1Lcorr'].array(library='np')
                data_y = np.take(data_y, np.where(data_y >= -6))
                mean_y = np.append(mean_y, data_y.astype('float64').mean())
                std_y = np.append(std_y, data_y.astype('float64').std())
            if data_type in ['DP1', 'DP2', 'DP']:
                label = r'{} runs'.format(data_type)
            elif data_type == 'neutron':
                label = r'Neutron runs'
            fluc_mean = max((max(mean_y)-np.mean(mean_y))/np.mean(mean_y), (min(mean_y)-np.mean(mean_y))/np.mean(mean_y))
            fluc_std = max((max(std_y)-np.mean(std_y))/np.mean(std_y), (min(std_y)-np.mean(std_y))/np.mean(std_y))
            print(r'Mean on log(S2/S1) {}: min {:.6f}, max {:.6f}, mean {:.6f}, variation {:.6f} %'.format(data_type,
                        min(mean_y), max(mean_y), np.mean(mean_y), fluc_mean*100))
            print(r'Std on log(S2/S1) {}: min {:.6f}, max {:.6f}, mean {:.6f}, variation {:.6f} %'.format(data_type,
                        min(std_y), max(std_y), np.mean(std_y), fluc_std*100))
            plt.errorbar(data_x, mean_y, yerr=std_y, fmt='.', capsize=6, elinewidth=1, label=label)
        plt.yticks(fontsize=font_size-2)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        plt.ylim(-1,1.6)
        f_indicate_date_ranges()
        plt.rcParams["legend.markerscale"] = 1.0
        plt.legend(fontsize=font_size-4, loc='lower left')
        if not plot_thesis_1col and not plot_thesis_2col:
            plt.title(r'Average of $\log_{10}(S2/S1)$ with standard deviation')
        plt.xlabel(r'Run date', fontsize=font_size+mag)
        plt.ylabel(r'S2maxLoS1Lcorr', fontsize=font_size+mag)
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()

    elif plot_name == 'S2maxLY_vs_RunDate':
        data_folder = ['DP', 'neutron']
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(alpha=0.5, zorder=0)
        print()
        for data_type in data_folder:
            data_x = np.array([], dtype=np.datetime64)
            mean_y = np.array([], dtype=np.float64)
            std_y = np.array([], dtype=np.float64)
            runs = next(os.walk('/mnt/raid/users/duylai/warios_' + data_type + '/'), (None, None, []))[2]
            summary = pd.read_csv('/mnt/raid/users/duylai/warios_' + data_type + '/satoshi-summary.txt', sep='\t')
            for run in runs:
                try:
                    run_nbr = int(run[4:8])
                    # run_nbr = int(run[8:12])
                except ValueError:
                    continue
                root_file = '/mnt/raid/users/duylai/Satoshi_' + data_type + '_runs/' + run
                # root_file = '/mnt/raid/users/duylai/warios_' + data_type + '/' + run
                tree = uproot.open(root_file + ':HebingTree')
                start_date = identify_dates(summary.query('SRunID  == ' + str(run_nbr))['StartRunID'].to_numpy()[0])
                end_date = identify_dates(summary.query('SRunID  == ' + str(run_nbr))['EndRunID'].to_numpy()[0])
                mid_date = np.datetime64(round(0.5 * (start_date.astype('int') + end_date.astype('int'))), 's')
                data_x = np.append(data_x, mid_date)
                data_y = tree['S2maxLY'].array(library='np')
                mean_y = np.append(mean_y, data_y.astype('float64').mean())
                std_y = np.append(std_y, data_y.astype('float64').std())
            if data_type in ['DP']:
                label = r'{} runs'.format(data_type)
            elif data_type == 'neutron':
                label = r'Neutron runs'
            print(r'Mean on S2max {}: min {:.2f}, max {:.2f}'.format(data_type, min(mean_y), max(mean_y)))
            plt.errorbar(data_x, mean_y, yerr=std_y, fmt='.', capsize=6, elinewidth=1, label=label)
        plt.yticks(fontsize=font_size-2)
        # ax.yaxis.set_major_locator(MultipleLocator(0.2))
        # plt.ylim(-1,1.6)
        f_indicate_date_ranges()
        plt.rcParams["legend.markerscale"] = 1.0
        plt.legend(fontsize=font_size-4, loc='lower left')
        plt.title(r'Average of S2maxL with standard deviation')
        plt.xlabel(r'Run date', fontsize=font_size+mag)
        plt.ylabel(r'S2maxL$', fontsize=font_size+mag)
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()

    elif plot_name == 'eff_detec_rec':
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(alpha=0.5, zorder=0)
        from quenching import *
        plot = uproot.open('/home/duylai/ardm/neutron/efficiency.root')
        energ = plot['h_detect;2'].to_numpy()[1]
        data = plot['h_detect;2'].to_numpy()[0]
        plt.plot(conversion_keV_pe(0.5*(energ[:-1]+energ[1:])), data
                , label=r'Detection')#, mean {:.4f}'.format(data.mean()))
        energ = plot['h_rec;2'].to_numpy()[1]
        data = plot['h_rec;2'].to_numpy()[0]
        plt.plot(conversion_keV_pe(0.5*(energ[:-1]+energ[1:])), data
                , label=r'Reconstruction')#, mean {:.4f}'.format(data.mean()))
        # StatsBox(False)
        plt.legend(fontsize=font_size + mag -2)
        plt.xlim(0,150)
        plt.ylim(0,1)
        plt.xlabel(r'S1L [p.e.]', fontsize=font_size+mag)
        plt.ylabel(r'Efficiency', fontsize=font_size+mag)
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()

    elif plot_name == 'e_depo':
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(alpha=0.5, zorder=0)
        from quenching import *
        plot = uproot.open('/home/duylai/ardm/neutron/efficiency.root')
        energ = plot['h_total;1'].to_numpy()[1]
        data = plot['h_total;1'].to_numpy()[0]
        plt.plot(conversion_keV_pe(0.5*(energ[:-1]+energ[1:])), data)
        # StatsBox(False)
        # plt.legend(fontsize=font_size + mag -2)
        plt.xlim(0,150)
        plt.yscale('log')
        plt.grid(which='major', ls='-')
        plt.grid(which='minor', ls='--', alpha=0.7)
        plt.ylim(1e1,1e5)
        plt.xlabel(r'S1L [p.e.]', fontsize=font_size+mag)
        plt.ylabel(r'Counts', fontsize=font_size+mag)
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()
    
    elif plot_name == 'eff_match':
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(alpha=0.5, zorder=0)
        energ = np.linspace(0,800,101)
        data = pd.read_csv('/home/duylai/ardm/eff_match.csv')['0']
        plt.plot(0.5*(energ[:-1]+energ[1:]), data)
        # StatsBox(False)
        plt.xlim(0,800)
        plt.ylim(0,1)
        plt.xlabel(r'S1L [p.e.]', fontsize=font_size+mag)
        plt.ylabel(r'Efficiency', fontsize=font_size+mag)
        ax.xaxis.set_major_locator(MultipleLocator(100.0))
        fig.tight_layout()
        if save_fig:
            plt.savefig(output_path + plot_name + file_type, dpi=300)
            plt.close()
        if show_fig:
            plt.show()