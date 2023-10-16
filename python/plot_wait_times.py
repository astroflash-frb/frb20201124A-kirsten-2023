import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from get_MJD import fit_gauss, multi_gauss
from matplotlib import rcParams
import os

rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['axes.labelsize'] = 10
rcParams['legend.fontsize'] = 8
rcParams['legend.handlelength'] = 1.0

db = '../dbs/burst_info.pickle'
data = pd.read_pickle(db)
data = data[(data.src == 'sfxc')]
ids = data.id.unique()

multicomp = []
waits = []

# multi-component bursts with detections at more than one dish
multicomp_multidish = [['B08-o8', 'B08-tr'],
                       ['B24-wb', 'B24-o8'],
                       ['B25-wb', 'B25-o8'],
                       ['B26-wb', 'B26-o8'],
                       ['B42-wb', 'B42-tr']]
# out of these we exclude the following in the plotting to avoid duplicates
exclude = ['B08-tr', 'B24-wb', 'B25-wb', 'B26-o8', 'B42-wb']


for ID in ids:
    if ID in exclude:
        continue
    ts = data[(data.id == ID)].toa_mjd_utc.values
    if len(ts) > 1:
        for wait in np.diff(ts):
            waits.append(wait * 86400)
        multicomp.append(ID)


n, bins, _ = plt.hist(waits, 50)

plt.figure(figsize=(5, 4))
plt.xlabel('Component separation [s]')
plt.savefig('../plots/wait_time_hist.png', dpi=300)

plt.clf()
waits = np.log10(waits)
# the below would create a log-log plot, but we only want log on the y-axis
# and so we just get the bins right first
n, bins, _ = plt.hist(waits, 12, log=True)
plt.clf()

plt.figure(figsize=(5, 4))
n, bins, _ = plt.hist(waits, bins)

# fit gaussian
amp = np.max(n)
sigma = bins.std()
x0 = np.median(bins)
params = np.array([[amp, x0, sigma]])
model = multi_gauss(params)

diffs = np.diff(bins) / 2
mbins = bins[:-1] + diffs
fit, fit_lm = fit_gauss(model, n, mbins)

plotbins = np.arange(bins.min(), bins.max(), 0.01)
plt.plot(plotbins, fit(plotbins))
plt.xlabel(r'$\mathrm{\log_{10}\left(Component\;separation / [second]\right)}$')
plt.ylabel('Count')
plt.savefig('../plots/wait_time_hist_log.png', dpi=300)
try:
    plt.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/wait_time_hist_log.png",
                    bbox_inches='tight',
                    facecolor='white', transparent=False, dpi=600)
except:
    print('Could not save to overleaf. Got only a local copy.')

print(fit)



# plot multi-component bursts detected at different dishes against one-another
'''
multicomp_multidish = [['B08-o8', 'B08-tr'],
                       ['B24-o8', 'B24-o8'],  # missing components in wb
                       ['B25-wb', 'B25-o8'],
                       ['B26-wb', 'B26-o8'],
                       ['B42-tr', 'B42-tr']] # missing components in o8
bursts = ['b08', 'b24', 'b25', 'b26', 'b42']
doubles = {}
for i, ids in enumerate(multicomp_multidish):
     for n, ID in enumerate(ids):
         print(n, ID)
         df = sfxc[(sfxc.id == ID) & (sfxc.src == 'sfxc')]
         ts = np.diff(df.toa_mjd_utc.values) * 86400 * 1e3
         if n == 0:
             bts = np.zeros((2, len(ts)), dtype=np.float)
         bts[n] = ts
     doubles[bursts[i]] = bts

for k in doubles.keys():
    values = doubles[k]
    plt.plot(values[0], values[1], 'x')
x = np.arange(2,12)
plt.plot(x,x)
'''
