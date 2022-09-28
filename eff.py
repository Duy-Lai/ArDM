# calculating the matching efficiency as a function of S1L by counting entries in .root files
# written by Duy Lai, 09.2022

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot, os

# names = ['matching']
names = ['whatever']

cut = 'preselection'

data_folder = ['neutron']
# data_folder = ['neutron', 'DP2']

for plot_name in names:
    if plot_name == 'matching':
        eff = np.array([], dtype=np.float64)
        slice = np.linspace(0,800,2)
        # eff = [0.5875045943067472, 0.5513633277837658, 0.54395966375306, 0.556787329787784, 0.5826299916384464,
        #         0.6159381451613746, 0.65090092772656, 0.6817857050686538, 0.7054933403862859, 0.7225513119900441]
        for i in np.arange(slice.shape[0]-1):
            satoshi = 0
            matched = 0
            for data in data_folder:
                pre_cut_path = '/mnt/raid/users/duylai/Satoshi_' + data + '_runs/'
                post_cut_path = '/mnt/raid/users/duylai/matched_' + data + '/'
                # print(next(os.walk(pre_cut_path), (None, None, []))[2])
                # print(next(os.walk(post_cut_path), (None, None, []))[2][:-1])
                for file in next(os.walk(pre_cut_path), (None, None, []))[2]:
                    try:
                        tree = uproot.open(pre_cut_path + file + ':HebingTree')
                    except ValueError:
                        continue
                    s1 = tree['S1LY'].array(library='np')
                    satoshi += np.count_nonzero((s1 >= slice[i]) & (s1 < slice[i+1]))
                for file in next(os.walk(post_cut_path), (None, None, []))[2][:-1]:
                    try:
                        tree = uproot.open(post_cut_path + file + ':HebingTree')
                    except ValueError:
                        continue
                    s1 = tree['S1LY'].array(library='np')
                    matched += np.count_nonzero((s1 >= slice[i]) & (s1 < slice[i+1]))
            try:
                print(matched)
                print(satoshi)
                eff = np.append(eff, matched / satoshi)
                # pd.Series(eff).to_csv('/home/duylai/ardm/neutron/eff_match_neutron.csv', index=False)
                print(r'Matching efficiency, S1L [{:.2f}, {:.2f}]: {}'.format(slice[i], slice[i+1], eff[i]*100))
            except ZeroDivisionError:
                continue

    elif plot_name == 'whatever':
        for data in data_folder:
            tree = uproot.open('/mnt/raid/users/duylai/' + cut + '_' + data + '/SuperWario.root:HebingTree')
            s2s1 = tree['S2maxLoS1Lcorr'].array(library='np')
            # s2dt = tree['S2maxDT'].array(library='np')[(s2s1 > -0.5) & (s2s1 <= 0.5)]
            s2dt = tree['S2maxDT'].array(library='np')[((s2s1 >= -3) & (s2s1 < -0.5)) | ((s2s1 > 0.5) & (s2s1 <= 5))]
            # plt.hist(s2dt)
            # plt.show()
            print(s2dt.shape[0])
            # s2dt_cut = s2dt[s2dt < 0.5]
            s2dt_cut = s2dt[s2dt >= 0.5]
            print(s2dt_cut.shape[0])
            print(s2dt_cut.shape[0]/s2dt.shape[0])
            # surv = np.count_nonzero(events >= 0.55)
            # print(surv/tree.member('fEntries')*100)
