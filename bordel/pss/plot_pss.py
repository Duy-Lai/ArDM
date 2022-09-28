import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator

types = ['dp_ER', 'neutron_NR']
bumps = ['upper', 'lower']

save_fig = 0

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
    return

fig, ax = plt.subplots()

plots = []
lgnd_text = []
for type in types:
    data = pd.read_csv('pss/pss_' + type + '.csv')
    dates = np.loadtxt('pss/run_dates_' + type + '.csv', dtype=np.datetime64)
    for bump in bumps:
        if type == 'dp_ER' and bump == 'upper':
            plots.append(plt.scatter(dates, np.array(data[bump], dtype=np.float64), c='b', marker='.', s=40))
            lgnd_text.append(r'Upper ER, $\log(S2/S1)\geq0.5$')
        elif type == 'dp_ER' and bump == 'lower':
            plots.append(plt.scatter(dates, np.array(data[bump], dtype=np.float64), c='g', marker='.', s=40))
            lgnd_text.append(r'Lower ER, $\log(S2/S1)<0.5$')
        elif type == 'neutron_NR' and bump == 'upper':
            plots.append(plt.scatter(dates, np.array(data[bump], dtype=np.float64), c='r', marker='x', s=40))
            lgnd_text.append(r'Upper NR, $\log(S2/S1)\geq0.2$')
        elif type == 'neutron_NR' and bump == 'lower':
            plots.append(plt.scatter(dates, np.array(data[bump], dtype=np.float64), c='orange', marker='x', s=40))
            lgnd_text.append(r'Lower NR, $\log(S2/S1)<0.2$')

f_indicate_date_ranges()
ax.grid(alpha=0.5, zorder=0)
plt.ylim(0,1.0)
ax.yaxis.set_major_locator(MultipleLocator(0.1))
lgnd = plt.legend(plots, lgnd_text)
lgnd.legendHandles[0]._sizes = [100]
lgnd.legendHandles[1]._sizes = [100]
lgnd.legendHandles[2]._sizes = [60]
lgnd.legendHandles[3]._sizes = [60]
plt.xlabel(r'Run date', fontsize=font_size+2)
plt.ylabel(r'Percentage of single scatters', fontsize=font_size+2)
plt.title(r'Single + multiple scattering',c='r',fontsize=font_size+2)
fig.tight_layout()
if save_fig == 1:
    plt.savefig('pss/pss_vs_runs.png', dpi=300)
else:
    plt.show()