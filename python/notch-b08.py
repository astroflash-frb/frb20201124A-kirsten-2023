import numpy as np
from presto import filterbank
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['axes.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['legend.handlelength'] = 1.0

#station = 'o8'
station = 'tr'

def get_spectra(filename, start, length, tscrunch, fscrunch):
    fil = filterbank.FilterbankFile(filname)
    tsamp = fil.header['tsamp'] # comes in s
    nchan = fil.header['nchans']
    dur = fil.nspec * fil.header['tsamp']
    f0 = fil.header['fch1']
    foff = fil.header['foff']
    fmax = f0 + np.abs(foff) / 2
    fmin = f0 + foff*nchan + np.abs(foff) / 2

    # back to second space
    secs_in_s = start / 1e3
    length_s = length / 1e3
    begbin = int(secs_in_s / tsamp)
    nspec = int(length_s / tsamp)
    spec = fil.get_spectra(begbin, nspec)
    offspec = fil.get_spectra(1000, 5260)

    spec.downsample(factor=tscrunch)
    offspec.downsample(factor=tscrunch)
    arr = spec.data
    offarr = offspec.data

    arr = (arr - offarr.mean(axis=1)[:,None]) / offarr.std(axis=1)[:, None]
    offarr = (offarr - offarr.mean(axis=1)[:,None]) / offarr.std(axis=1)[:, None]
    if fscrunch > 1:
        arr = arr.reshape(-1, fscrunch, arr.shape[1]).mean(axis=1)
        offarr = offarr.reshape(-1, fscrunch, offarr.shape[1]).mean(axis=1)
    # getting the time series as well
    ts = arr.mean(axis=0)
    off_ts = offarr.mean(0)
    ts -= off_ts.mean()
    ts /= off_ts.std()
    return arr, ts, fmin, fmax


if station == 'o8':
    filname = '../data/r67028_o8_no0020_sfxc_cutout_ds.fil'
    rect1 = patches.Rectangle((7.5, 35), 5.5, -45, linewidth=2, edgecolor='r', facecolor='none')
    rect2 = patches.Rectangle((7.5, 1360), 5.5, 132, linewidth=2, edgecolor='r', facecolor='none')
elif station == 'tr':
    filname = '../data/r67l01_tr_no0029_sfxc_cutout_ds.fil'
    rect1 = patches.Rectangle((7.5, 35), 5.5, -45, linewidth=2, edgecolor='r', facecolor='none')
    rect2 = patches.Rectangle((7.5, 1340), 5.5, 260, linewidth=2, edgecolor='r', facecolor='none')

else:
    print(f'Unknown station {station}. Only o8 and tr supported.')
    quit(1)

# all times in ms
zoom_start = 115
zoom_length = 6
full_start = 108
full_length = 31
z_tscrunch = 4
f_tscrunch = 16
z_fscrunch = 4
f_fscrunch = 4

zoom, zoom_ts, fmin, fmax = get_spectra(filname, start=zoom_start, length=zoom_length,
                               tscrunch=z_tscrunch, fscrunch=z_fscrunch)
full, full_ts, fmin, fmax = get_spectra(filname, start=full_start, length=full_length,
                               tscrunch=f_tscrunch, fscrunch=f_fscrunch)
#zoom_ts = zoom.mean(axis=0)
#full_ts = full.mean(axis=0)

fig = plt.figure(figsize=(4, 4))
rows = 2
cols = 2
gs = gridspec.GridSpec(rows, cols, wspace=0.15, hspace=0,
                       height_ratios=[1, 3],
                       width_ratios=[1, 1])

length = full_length
profile_full = full_ts
arr = full
x_axis = np.linspace(0, length, len(profile_full))
prof = plt.subplot(gs[0])
prof.plot(x_axis, profile_full)
prof.set_ylabel('S/N')
prof.tick_params(labelbottom=False)
#prof.set_xlim(pltrange)
vmin = arr.std() * -1
vmax = arr.std() * 6
dyn = plt.subplot(gs[2], sharex=prof)
dyn.imshow(arr, aspect='auto', interpolation='None',
           vmin=vmin, vmax=vmax, extent=(x_axis[0], x_axis[-1], fmin, fmax))
#dyn.set_xlim(prof.get_xlim())
xmax, xmin = dyn.get_xlim()
middle = (xmin + xmax) // 2
nticks = 7
xtick_lables = np.linspace(-length//2, length//2, nticks, dtype=int).astype(str)
xticks = np.linspace(middle-length//2, middle+length//2, nticks)
dyn.set_xticks(xticks)
dyn.set_xticklabels(xtick_lables)
dyn.set_ylabel('Frequency [MHz]')
dyn.set_xlabel('Time [ms]')

length2 = zoom_length
profile_full2 = zoom_ts
arr2 = zoom
x_axis2 = np.linspace(0, length2, len(profile_full2))
prof2 = plt.subplot(gs[1])
prof2.plot(x_axis2, profile_full2)
prof2.tick_params(labelbottom=False)

vmin = arr2.std() * -1
vmax = arr2.std() * 6
dyn2 = plt.subplot(gs[3], sharex=prof2)
dyn2.imshow(arr2, aspect='auto', interpolation='None',
           vmin=vmin, vmax=vmax, extent=(x_axis2[0], x_axis2[-1], fmin, fmax))

xmax, xmin = dyn2.get_xlim()
middle = (xmin + xmax) // 2
nticks = 7
xtick_lables = np.linspace(-length2//2, length2//2, nticks, dtype=int).astype(str)
xticks = np.linspace(middle-length2//2, middle+length2//2, nticks)
dyn2.set_xticks(xticks)
dyn2.set_xticklabels(xtick_lables)
dyn2.set_xlabel('Time [ms]')
dyn2.tick_params(labelleft=False)

# Add the patch to the Axes
prof.add_patch(rect1)
dyn.add_patch(rect2)

plt.savefig(f'../plots/notch-b08-{station}.png', bbox_inches='tight',
            facecolor='white', transparent=False, dpi=600)
try:
    plt.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/notch-b08-{station}.png",
                bbox_inches='tight',
                facecolor='white', transparent=False, dpi=600)
except:
    print('Could not save plot to overleaf repo. Got a local copy only.')
