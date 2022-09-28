# used for the quantification of source events, requires CM_results.csv produced by Fit2D.cpp
# written by Duy Lai, 09.2022

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import (MultipleLocator, ScalarFormatter, LogLocator)
import os, sys, uproot, json

sys.path.insert(0, '/home/duylai/ardm/')
from __Var__ import *
from quenching import *

################################ KNOWN ACTIVITY ################################

satoshi = pd.read_csv('/mnt/raid/users/duylai/Satoshi_neutron_runs/satoshi-summary.txt', delimiter='\t')
All_runs = pd.read_csv('All_runs.csv')

# query first and last RunIDs
for i in np.arange(satoshi.shape[0]):
    run = satoshi.iloc[i,:]
    if run.Type in ['DP2', 'unk']:
        start_time = np.datetime64(All_runs.query('RunID == {}'.format(run.StartRunID)).Date.to_numpy()[0], 's')
        end_time = np.datetime64(All_runs.query('RunID == {}'.format(run.EndRunID)).Date.to_numpy()[0], 's')
        # total_runtime += np.float64(end_time - start_time)
        # print('Run {}: from {} to {} -> {:.1f} s'.format(run.SRunID, start_time, end_time, np.float64(end_time-start_time)))
        # nruns += (run.EndRunID-run.StartRunID)
        # nruns += run['#Runs']
        try:
            RunID_min = np.minimum(run.StartRunID, RunID_min)
            RunID_max = np.maximum(run.EndRunID, RunID_max)
        except NameError:
            RunID_min = run.StartRunID
            RunID_max = run.EndRunID

def cf252_activity(date):
    neutrons_per_SF = 3.7692
    SF_branching_ratio = 0.03092
    activity_1816 = 1e4
    half_life = 2.645
    decay_rate = np.log(2) / (half_life * 365.2421987)
    return neutrons_per_SF * SF_branching_ratio * activity_1816 * np.exp(-decay_rate * date)

# obtaining the time/date of the first run and the last
start_time = np.datetime64(All_runs.query('RunID == {}'.format(RunID_min)).Date.to_numpy()[0], 's')
end_time = np.datetime64(All_runs.query('RunID == {}'.format(RunID_max)).Date.to_numpy()[0], 's')

# calculating the corresponding activities on the first date and the last, grouped into a list of 2 elements
activity_th = cf252_activity(np.float64([start_time, end_time]-np.datetime64('2016-08-01', 's'))/3600/24)

# activity is taken as the mean between the first and the last, error is their difference divided by 2
activity_known = Var(np.mean(activity_th), np.mean(activity_th * [1,-1]))

print('{}:\n{}'.format('²⁵²Cf activity', activity_known))

################################ RUN TIME ################################

v3 = next(os.walk('satoshiV3/'), (None, None, []))[2]

nruns = 0
total_runtime = 0

for v in v3:
    runs = json.load(open('satoshiV3/' + v))['runIDs']
    nruns += len(runs)
    for r in runs:
        total_runtime += All_runs.query('RunID == {}'.format(r)).DAQTimeLapse.to_numpy()[0]

# print(total_runtime)
# print(nruns)

run_time = Var(total_runtime / nruns, 1.0) * 3688
# run_time = Var(115791.0, 3688.0)

print('Total runtime:\n{}'.format(run_time))

# activity * total run time
N_source_cf = VarProd([activity_known, run_time])

print('Expected number of emitted ²⁵²Cf events:\n{}'.format(N_source_cf))

# acceptance efficiency = 0.07
eff_geom = Var(0.07, 0.0)   # preliminary
N_source_interac = VarProd([N_source_cf, eff_geom])

print('Expected number of interacted ²⁵²Cf events:\n{}'.format(N_source_interac))

################################ MONTE-CARLO & EFFICIENCIES ################################

MC_file = uproot.open('neutron/efficiency.root')
range = np.linspace(0,100,11)

val_hist = MC_file['h_total;1'].to_numpy()
pe_range = conversion_keV_pe(0.5 * (val_hist[1][:-1] + val_hist[1][1:]))
# multiplying fraction of MC events below 100 p.e., N_source_interac has 1 element
N_source_interac *= np.sum(val_hist[0][pe_range <= 100]) / np.sum(val_hist[0])

# print(val_hist[0].sum())
# print(MC_file['h_total;1'].member('fEntries'))
# print(N_source_interac)

e_depo = np.array([])
for i in np.arange(range.shape[0] - 1):
    e_depo = np.append(e_depo, np.sum(val_hist[0][(pe_range >= range[i]) & (pe_range < range[i+1])]))

# scaling the MC spectrum to the "real" number of events, N_source_interac now has 10 elements
N_source_interac *= (e_depo / e_depo.sum())

# print(val_hist[0][pe_range<100].sum())
# print(e_depo.sum())

# print(N_source_interac.sum_val()*(1/e_depo.sum()))

# Detection efficiency
val_hist = MC_file['h_detect;2'].to_numpy()
pe_range = conversion_keV_pe(0.5 * (val_hist[1][:-1] + val_hist[1][1:]))
mean_eff = np.array([])
std_eff = np.array([])
for i in np.arange(range.shape[0] - 1):
    mean_eff = np.append(mean_eff, np.mean(val_hist[0][(pe_range >= range[i]) & (pe_range < range[i+1])]))
    std_eff = np.append(std_eff, np.std(val_hist[0][(pe_range >= range[i]) & (pe_range < range[i+1])]))

eff_detect = Var(mean_eff, std_eff)

# plt.plot(np.linspace(0,100,10), mean_eff)
# plt.grid()
# plt.show()

# Reconstruction efficiency
val_hist = MC_file['h_rec;2'].to_numpy()
pe_range = conversion_keV_pe(0.5 * (val_hist[1][:-1] + val_hist[1][1:]))
mean_eff = np.array([])
std_eff = np.array([])
for i in np.arange(range.shape[0] - 1):
    mean_eff = np.append(mean_eff, np.mean(val_hist[0][(pe_range >= range[i]) & (pe_range < range[i+1])]))
    std_eff = np.append(std_eff, np.std(val_hist[0][(pe_range >= range[i]) & (pe_range < range[i+1])]))

eff_recon = Var(mean_eff, std_eff)

# Matching efficiency, eff_match_neutron.csv produced by eff.py
eff_match = Var(pd.read_csv('eff_match_neutron.csv')['0'].to_numpy(), np.zeros(10))

# eff_match = Var([0.5875045943067472, 0.5513633277837658, 0.54395966375306, 0.556787329787784, 0.5826299916384464,
#         0.6159381451613746, 0.65090092772656, 0.6817857050686538, 0.7054933403862859, 0.7225513119900441], np.zeros(10))

# MC spectrum of matched events
N_source_prem = VarProd([N_source_interac, eff_detect, eff_recon, eff_match])

print('{}:\n{}'.format('Number of pre-matched ²⁵²Cf events', N_source_prem))

################################ SOURCE EVENTS' EFFICIENCY ################################

# Number of events in reality, resultes of Fit2D.cpp
cm = pd.read_csv('neutron/CM_results.csv')

N_ER_DP = Var(cm.N_ER_DP, cm.N_ER_DP_err)
N_NR_DP = Var(cm.N_NR_DP, cm.N_NR_DP_err)
N_DP = VarSum([N_ER_DP, N_NR_DP])

N_ER_Neutron = Var(cm.N_ER_Neutron, cm.N_ER_Neutron_err)
N_NR_Neutron = Var(cm.N_NR_Neutron, cm.N_NR_Neutron_err)
N_Neutron = VarSum([N_ER_Neutron, N_NR_Neutron])

N_instr_DP = N_NR_DP
frac_instr_DP = VarProd([N_instr_DP, N_DP], [1])

# instrumental background
N_instr_Neutron = VarProd([N_instr_DP, N_ER_Neutron, N_ER_DP], [2])
frac_instr_Neutron = VarProd([N_instr_Neutron, N_Neutron], [1])

# number of source events in reality
N_source = VarSum([N_NR_Neutron, N_instr_Neutron], [1])

print('{}:\n{}'.format('Number of Fit-2D ²⁵²Cf events', N_source))

# N_source_prem_sum = N_source_prem.sum_val()

# print('{}:\n{}'.format('Sum of pre-matched ²⁵²Cf events', N_source_prem_sum))

# Pre-selection and cuts of interest efficiency for source events
eff_cuts = VarProd([N_source, N_source_prem], [1])

print('Cut efficiencies: \n{}'.format(eff_cuts))

# print('Sum of the cut efficiencies: {}'.format(eff_cuts.sum_val()))

################################ PLOTTING ################################

def make_error_boxes(ax, y, xdata, xerror, facecolor='b',
                     edgecolor='none', alpha=0.2):
    ydata = y.val
    yerror = y.err
    # Loop over fit2D points; create box from errors at each point
    errorboxes = [Rectangle((x - xe, y - ye), xe*2, ye*2)
                  for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T)]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)
    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror, capsize=2
                        #   , fmt='o', color='k', ecolor=facecolor)
                          , fmt='o', color='k', ecolor='k')
    return artists

mag = 2
# mag = 6

# Plotting parameters
font_size = 14
plt.rcParams["axes.labelsize"] = font_size +2
plt.rcParams["legend.framealpha"] = 1.0
# plt.rcParams["axes.titlesize"] = title_size
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = font_size - 4
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams["legend.markerscale"] = 1
plt.rc('xtick', labelsize=font_size + mag)
plt.rc('ytick', labelsize=font_size + mag)

# fig, ax = plt.subplots(figsize=(7,5))
fig, ax = plt.subplots(figsize=(8,5))
make_error_boxes(ax, eff_cuts, (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy())
# make_error_boxes(ax, N_source_interac, (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy())
# make_error_boxes(ax, N_source_prem, (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy(),
    # facecolor='r')
# make_error_boxes(ax, VarProd([N_source_cf, eff_geom, eff_detect]), (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy(),
#     facecolor='r')
# make_error_boxes(ax, VarProd([N_source_cf, eff_geom, eff_detect, eff_match]), (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy(),
#     facecolor='g')
# make_error_boxes(ax, VarProd([N_source_cf, eff_geom, eff_detect, eff_match, eff_match]), (0.5*(cm.S1L_low + cm.S1L_high)).to_numpy(), (0.5*(cm.S1L_high - cm.S1L_low)).to_numpy(),
#     facecolor='b')
# plt.fill_between(x, y + yerr, y - yerr, color='b', alpha=0.2)
# plt.legend([r'$N_{\rm source}^{\rm interacted}$',r'$N_{\rm source}^{\rm matched}$'],
        #  fontsize = font_size + mag +2)
plt.grid(which='major', ls='-')
plt.grid(which='minor', ls='--', alpha=0.7)
plt.xlabel('S1L [p.e.]', fontsize = font_size + mag)
plt.ylabel('Cut efficiency', fontsize = font_size + mag)
# plt.ylabel('Expected number of events', fontsize = font_size + mag)
plt.xlim(0,100)
ax.xaxis.set_major_locator(MultipleLocator(10.0))
# plt.yscale('log')
# plt.ylim(1e3,1e7)
plt.ylim(0,0.16)
# ax.yaxis.set_major_locator(MultipleLocator(200.0))
# ax.yaxis.set_major_locator(LogLocator(numticks=15))
# formatter = ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1,1))
# ax.yaxis.set_major_formatter(formatter)
fig.tight_layout()
plt.show()
# plt.savefig('/home/duylai/ardm//thesis_plots/eff_cuts.svg', dpi=300)