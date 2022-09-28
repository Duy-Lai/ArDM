import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import h5py, os, json

recoil_type = 'NR'
save_fig = 0

if save_fig == 0:
    show_fig = 1
else:
    show_fig = 0

dp_file = pd.read_csv('figures_dp_' + recoil_type + '/runs/entries_bumps.csv')
neutron_file = pd.read_csv('figures_neutron_' + recoil_type + '/runs/entries_bumps.csv')

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

def identify_date(source):
    if source == 'dp':
        file = dp_file
    elif source == 'neutron':
        file = neutron_file
    else:
        raise ValueError('Invalid recoil type')
    run_IDs = np.array([], dtype = int)
    for i in file.run:
        for s in satoshi:
            if int(s[5:8]) == i and s[9:12] != 'unk' and ((s[13:17] == 'none' and source == 'dp')
                    or (s[13:20] == 'neutron' and source == 'neutron')):
                run_list = json.load(open('satoshiV3/' + s))['runlist']
        run_IDs = np.append(run_IDs, np.where(fh5_runIDs == run_list[round(len(run_list) / 2)]))
    return fh5_dates[run_IDs], len(run_list)

run_dates_dp, nruns_dp = identify_date('dp')
run_dates_neutron, nruns_neutron = identify_date('neutron')

# np.savetxt('run_dates_nr.csv', run_dates_nr, fmt='%s')

#####################################################################################################

# Plotting parameters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size + 2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size - 4
plt.rcParams["legend.framealpha"] = 1.0
plt.rcParams["legend.markerscale"] = 3.5
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)

# Indicating specific dates
def f_indicate_date_ranges():
    # locator = mdates.AutoDateLocator()
    locator = mdates.DayLocator(bymonthday=(1, 6, 11, 16, 21, 26))
    formatter = mdates.ConciseDateFormatter(locator, show_offset=False)
    plt.xlim(np.datetime64('2019-06-11'), np.datetime64('2019-09-01'))
    date_proxy = np.arange(np.datetime64('2019-05-12T14:35:49'), np.datetime64('2019-08-26T12:51:30'),
                np.timedelta64(30, 's'))
    plt.xticks(rotation=50)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

fig, ax = plt.subplots()
p1 = plt.scatter(run_dates_neutron,neutron_file.upper/nruns_neutron,marker='x',s=40,color='r')
p2 = plt.scatter(run_dates_dp,dp_file.upper/nruns_dp,s=7,color='b')
p3 = plt.scatter(run_dates_neutron,neutron_file.lower/nruns_neutron,marker='x',s=40,color='orange')
p4 = plt.scatter(run_dates_dp,dp_file.lower/nruns_dp,s=7,color='g')

plt.xlabel(r'Run date',fontsize=font_size+2)
plt.ylabel(r'Events per 1 run',fontsize=font_size+2)
# plt.ylabel(r'Fraction of events',fontsize=font_size+2)

f_indicate_date_ranges()
# plt.ylim(0,1)
# ax.yaxis.set_major_locator(MultipleLocator(0.1))

if recoil_type == 'NR':
    plt.title(r'NR events (RoI)')
elif recoil_type == 'ER':
    plt.title(r'ER events (Inverse RoI)')

lgnd = plt.legend([p2,p4,p1,p3], [r'DP, $\log(S2/S1)\geq0.5$', r'DP, $\log(S2/S1)<0.5$',
            r'Neutron, $\log(S2/S1)\geq0.2$', r'Neutron, $\log(S2/S1)<0.2$'], loc=0)
lgnd.legendHandles[0]._sizes = [40]
lgnd.legendHandles[1]._sizes = [40]
lgnd.legendHandles[2]._sizes = [40]
lgnd.legendHandles[3]._sizes = [40]

plt.grid(alpha=0.5, zorder=0)
fig.tight_layout()

if save_fig == 1:
    plt.savefig('bump_ratios/scattering_' + scatter + '.png', dpi=300)
else:
    plt.show()