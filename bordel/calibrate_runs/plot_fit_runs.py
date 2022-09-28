import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import h5py, os, json, math

y_col = 5
# possible choices for y_col: {3, 5, 7, 8}
bump_nbr = 1
as_fct_of ='runs'
save_fig = 0

if save_fig == 0:
    show_fig = 1
else:
    show_fig = 0

if y_col not in (3,5,7,8):
    raise ValueError('Choose (2 or 4 or 6 or 7) for y_col')

if bump_nbr not in (1,2):
    raise ValueError('Choose (1 or 2) for bump_nbr')

#####################################################################################################

# Read data from the fits
if bump_nbr == 1:
    runs_single_nr = pd.read_csv('calibrate_' + as_fct_of + '/neutron_NR/single_upper/fit_data.csv')
    runs_multiple_nr = pd.read_csv('calibrate_' + as_fct_of + '/neutron_NR/multi_upper/fit_data.csv')
    runs_single_er = pd.read_csv('calibrate_' + as_fct_of + '/dp_ER/single_upper/fit_data.csv')
    runs_multiple_er = pd.read_csv('calibrate_' + as_fct_of + '/dp_ER/multi_upper/fit_data.csv')
elif bump_nbr == 2:
    runs_single_nr = pd.read_csv('calibrate_' + as_fct_of + '/neutron_NR/single_lower/fit_data.csv')
    runs_multiple_nr = pd.read_csv('calibrate_' + as_fct_of + '/neutron_NR/multi_lower/fit_data.csv')
    runs_single_er = pd.read_csv('calibrate_' + as_fct_of + '/dp_ER/single_lower/fit_data.csv')
    runs_multiple_er = pd.read_csv('calibrate_' + as_fct_of + '/dp_ER/multi_lower/fit_data.csv')

nruns_nr = [30, 52, 53, 54, 56, 60, 63, 64, 65, 66, 84, 85, 100, 115]
nruns_er = []
for er in np.arange(125):
    if er not in nruns_nr:
        nruns_er = np.append(nruns_er, er)

# Obtain all the dates and their corresponding IDs
file = 'All_runs.h5'
fh5 = h5py.File(file, 'r')
fh5_runIDs = np.array([])
fh5_runIDs = np.append(fh5_runIDs, fh5['runID'][:])
fh5_dates = np.array([], dtype = np.datetime64)
for t in fh5['DAQStartTime'][:]:
    t = t.decode('utf-8')
    dt = np.array([t[6:10] + '-' + t[0:2] + '-' + t[3:5] + 'T' + t[12:20]], dtype = np.datetime64)
    fh5_dates = np.append(fh5_dates, dt)

# Retrieve the list of satoshi-runs
satoshi = next(os.walk('satoshiV3'), (None, None, []))[2]

# Match the run IDs with actual dates
def identify_dates(scattering, recoil):
    if scattering not in ('single', 'multi') or recoil not in ('NR', 'ER'):
        raise ValueError('Invalid scatering and/or recoil type')
        
    if scattering == 'single' and recoil == 'NR':
        nruns = nruns_nr
        csv_file = runs_single_nr
    elif scattering == 'single' and recoil == 'ER':
        nruns = nruns_er
        csv_file = runs_single_er
    elif scattering == 'multi' and recoil == 'NR':
        nruns = nruns_nr
        csv_file = runs_multiple_nr
    elif scattering == 'multi' and recoil == 'ER':
        nruns = nruns_er
        csv_file = runs_multiple_er

    run_IDs = np.array([], dtype = int)
    start_IDs = np.array([], dtype = int)
    end_IDs = np.array([], dtype = int)

    for i in np.arange(csv_file.start.size):
        run_list = np.array([], dtype = int)
        for j in nruns[np.where(nruns == csv_file.iloc[i,0])[0][0]:(np.where(nruns == csv_file.iloc[i,1])[0]+1)[0]]:
            for s in satoshi:
                if int(s[5:8]) == j and s[9:12] != 'unk' and ((s[13:20] == 'neutron' and recoil == 'NR')
                        or (s[13:17] == 'none' and recoil == 'ER')):
                    run_list = np.append(run_list, json.load(open('satoshiV3/' + s))['runlist'])
        run_IDs = np.append(run_IDs, np.where(fh5_runIDs == run_list[round(len(run_list) / 2)]))
        start_IDs = np.append(start_IDs, np.where(fh5_runIDs == run_list[0]))
        end_IDs = np.append(end_IDs, np.where(fh5_runIDs == run_list[-1]))

    return fh5_dates[run_IDs], fh5_dates[start_IDs], fh5_dates[end_IDs]

run_single_nr_dates, start_single_nr_dates, end_single_nr_dates = identify_dates('single', 'NR')
run_single_er_dates, start_single_er_dates, end_single_er_dates = identify_dates('single', 'ER')
run_multi_nr_dates, start_multi_nr_dates, end_multi_nr_dates = identify_dates('multi', 'NR')
run_multi_er_dates, start_multi_er_dates, end_multi_er_dates = identify_dates('multi', 'ER')

#####################################################################################################

# Plotting parameters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size + 2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams["legend.markerscale"] = 3.5
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)

# Indicating specific dates
def f_indicate_date_ranges():
    # locator = mdates.AutoDateLocator()
    locator = mdates.DayLocator(bymonthday=(1, 5, 9, 13, 17, 21, 25, 29))
    formatter = mdates.ConciseDateFormatter(locator, show_offset=False)
    plt.xlim(np.datetime64('2019-06-17'), np.datetime64('2019-08-25'))
    date_proxy = np.arange(np.datetime64('2019-05-12T14:35:49'), np.datetime64('2019-08-26T12:51:30'),
                np.timedelta64(30, 's'))
    plt.xticks(rotation=50)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    return

# Fractions
fig, ax = plt.subplots()

if y_col == 8:
    p1, = plt.plot(run_single_nr_dates, runs_single_nr.iloc[:,y_col]*180/math.pi, marker='x', markersize=8,
                linewidth=1.5, color='b')
    p2, = plt.plot(end_multi_nr_dates, runs_multiple_nr.iloc[:,y_col]*180/math.pi, marker='x', markersize=8,
                linewidth=1.5, color='r')
    p3, = plt.plot(run_single_er_dates, runs_single_er.iloc[:,y_col]*180/math.pi, marker='o', markersize=8,
                linewidth=1.5, color='g')
    p4, = plt.plot(end_multi_er_dates, runs_multiple_er.iloc[:,y_col]*180/math.pi, marker='o', markersize=8,
                linewidth=1.5, color='orange')
else:
    p1, = plt.plot(run_single_nr_dates, runs_single_nr.iloc[:,y_col], marker='x', markersize=8,
                linewidth=1.5, color='b')
    p2, = plt.plot(end_multi_nr_dates, runs_multiple_nr.iloc[:,y_col], marker='x', markersize=8,
                linewidth=1.5, color='r')
    p3, = plt.plot(run_single_er_dates, runs_single_er.iloc[:,y_col], marker='*', markersize=8,
                linewidth=1.5, color='g')
    p4, = plt.plot(run_multi_er_dates, runs_multiple_er.iloc[:,y_col], marker='*', markersize=8,
                linewidth=1.5, color='orange')
if y_col in (3,5):
    f1 = plt.fill_between(run_single_nr_dates, runs_single_nr.iloc[:,y_col]+runs_single_nr.iloc[:,y_col+1],
                runs_single_nr.iloc[:,y_col]-runs_single_nr.iloc[:,y_col+1], color='b', alpha=0.25, hatch='\\')
    f2 = plt.fill_between(end_multi_nr_dates, runs_multiple_nr.iloc[:,y_col]+runs_multiple_nr.iloc[:,y_col+1],
                runs_multiple_nr.iloc[:,y_col]-runs_multiple_nr.iloc[:,y_col+1], color='r', alpha=0.25, hatch='\\')
    f3 = plt.fill_between(run_single_er_dates, runs_single_er.iloc[:,y_col]+runs_single_er.iloc[:,y_col+1],
                runs_single_er.iloc[:,y_col]-runs_single_er.iloc[:,y_col+1], color='g', alpha=0.25, hatch='//')
    f4 = plt.fill_between(run_multi_er_dates, runs_multiple_er.iloc[:,y_col]+runs_multiple_er.iloc[:,y_col+1],
                runs_multiple_er.iloc[:,y_col]-runs_multiple_er.iloc[:,y_col+1], color='orange', alpha=0.25, hatch='//')
## plt.errorbar(run_single_nr_dates, runs_single_nr.iloc[:,y_col], runs_single_nr.sigma_F90, (end_single_nr_dates-start_single_nr_dates)/2.0, marker='o',
##             markersize=4, capsize=8.0, alpha=0.7, linewidth=2.5, label='Single scattering')
## plt.errorbar(run_single_nr_dates, runs_multiple_nr.iloc[:,y_col], runs_multiple_nr.sigma_F90, (end_single_nr_dates-start_single_nr_dates)/2.0, marker='o',
##              markersize=4, capsize=8.0, alpha=0.7, linestyle='None', label='Multiple scattering')
plt.xlabel(r'Run date', fontsize=font_size+2)
f_indicate_date_ranges()

if bump_nbr == 1:
    plt.title('Upper bump', fontsize=font_size+2, color='r')
elif bump_nbr == 2:
    plt.title('Lower bump', fontsize=font_size+2, color='r')

lgnd_loc = 0

if y_col == 3:
    plt.ylabel(r'F90', fontsize=font_size+2)
    plt.ylim(0.0, 1.0)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    lgnd_loc = 6
elif y_col == 5:
    plt.ylabel(r'$\log_{10}\left(\frac{S2}{S1}\right)$', fontsize=font_size+2)
    plt.ylim(-3.0, 4.0)
    #ax.yaxis.set_major_locator(MultipleLocator(0.5))
elif y_col == 7:
    plt.ylabel(r'Correlation coefficient', fontsize=font_size+2)
    plt.ylim(-1.0, 1.0)
elif y_col == 8:
    plt.ylabel(r'$\theta_{Gauss.}$ [deg.]', fontsize=font_size+2)
    plt.ylim(-10.0, 10.0)

if y_col in (3,5):
    lgnd = plt.legend([(p1,f1), (p2,f2), (p3, f3), (p4, f4)], [r'NR single', r'NR multiple', r'ER single', r'ER multiple'],
                framealpha=0.5, loc=lgnd_loc)
elif y_col in (7,8):
    lgnd = plt.legend([p1, p2], [r'Single scattering', r'Multiple scattering'])
lgnd.legendHandles[0].set_markersize(8)
lgnd.legendHandles[0].set_linewidth(1.5)
lgnd.legendHandles[1].set_markersize(8)
lgnd.legendHandles[1].set_linewidth(1.5)
lgnd.legendHandles[2].set_markersize(8)
lgnd.legendHandles[2].set_linewidth(1.5)
lgnd.legendHandles[3].set_markersize(8)
lgnd.legendHandles[3].set_linewidth(1.5)

ax.grid(alpha=0.5, zorder=0)
fig.tight_layout()

file_name = 'calibrate_' + as_fct_of + '/fit_' + runs_single_nr.columns[y_col] + '_vs_' + as_fct_of
if bump_nbr == 1:
    file_name += '_upper'
elif bump_nbr == 2:
    file_name += '_lower'
file_name += '.png'

if save_fig:
    plt.savefig(file_name, dpi=300)
if show_fig:
    plt.show()

plt.close()
del fig, ax
