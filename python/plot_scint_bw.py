import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

db = '../dbs/burst_info.pickle'
data = pd.read_pickle(db)
data = data[(data.src == 'sfxc')]
ids = data.id.unique()

multicomp = []
waits = []

for ID in ids:
    ts = data[(data.id == ID)].toa_mjd_utc.values
    if len(ts) > 1:
        for wait in np.diff(ts):
            waits.append(wait * 86400 * 1e3)
        multicomp.append(ID)

n, bins, _ = plt.hist(waits, 50)

plt.xlabel('Component separation [ms]')
plt.savefig('../plots/wait_time_hist.png', dpi=300)
