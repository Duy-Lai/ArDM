import numpy as np
import pandas as pd
import uproot, os

types = ['dp_ER', 'neutron_NR']
bumps = ['', '1bump/', '2bump/']

for type in types:
    run_dates = np.loadtxt('/home/duylai/ardm/pss/run_dates_' + type + '.csv', dtype=np.datetime64)
    data = pd.DataFrame(index=np.arange(len(run_dates)))
    for bump in bumps:
        one_column = np.array([], dtype=np.float64)
        file_list = next(os.walk('/mnt/raid/users/duylai/warios_' + type + '/' + bump), (None, None, []))[2]
        for file in file_list:
            try:
                int(file[8:12])
                tree = uproot.open('/mnt/raid/users/duylai/warios_' + type + '/' + bump + file + ':HebingTree')
                S2maxFrac = np.array(tree['S2maxFracS2L'], dtype=np.float64)
                one_column = np.append(one_column, np.count_nonzero(S2maxFrac==1.0)/S2maxFrac.shape[0])
                # data = data.append(pd.DataFrame([[pss]], columns=[bump]), ignore_index=False)
            except ValueError:
                pass
        data = data.join(pd.DataFrame({bump: one_column}))
    data.rename(columns={'':'both', '1bump/':'upper', '2bump/':'lower'}, inplace=True)
    data.to_csv('/home/duylai/ardm/pss/pss_' + type + '.csv', index=False)
