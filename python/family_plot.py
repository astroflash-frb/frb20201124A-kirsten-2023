#!/usr/bin/env python3
# coding: utf-8
from load_file_kenzie import load_archive
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LinearLocator
import matplotlib as mpl
import os,copy
from tqdm import tqdm
import numpy.ma as ma
import pandas as pd
from station_colors import get_station_colors


def subband_edges(bandwidth, n_chan, tele):
    if tele == "st":
        subband_edge = n_chan // 1
        ab = np.array([0,1,subband_edge-2,subband_edge-1])
        print("bandwidth is:", bandwidth, ". Number of channels:", n_chan, ". First subband array: ", ab)
        print("Subband edges: ", subband_edge)

        mask_arr = ab

    elif bandwidth == 512 and tele == 'o8':
        subbands = 16
        mask_arr = subband_mask(bandwidth, n_chan, subbands)

    elif bandwidth == 256 and tele == 'o8':
        subbands = 16
        mask_arr = subband_mask(bandwidth, n_chan, subbands)

    else:
        subbands=8
        mask_arr = subband_mask(bandwidth, n_chan, subbands)
    return mask_arr


def subband_mask(bandwidth, n_chan, subbands):

    # ab is starting array [0,1,126,127]
    print(n_chan)
    subband_edge = n_chan // subbands

    #if statement so to make sure then bandpass also gets smaller when downsampling
    ab = np.array([0,1,subband_edge-2,subband_edge-1])

    if n_chan == 512:
        ab = np.array([0,subband_edge-1])

    if n_chan <= 256:
        ab = np.array([0])

    print("bandwidth is:", bandwidth, ". Number of channels:", n_chan, ". First subband array: ", ab)
    print("Subband edges: ", subband_edge)

    mask = []

    for i in range(subbands):
        a = ab+(i*subband_edge)
        mask.append(a)

    mask_arr = np.asarray(mask).ravel()

    return mask_arr


def normalise(ds, t_cent, bw, n_chan, tele, on_mask=None, t_sig=None):
    """
    Calibrate the dynamic spectrum for the bandpass

    Per frequency channel it subtracts the off burst mean and divides by the off burst std
    """

    # on_mask is True where things are masked
    ds_off = np.ma.masked_array(ds, mask=on_mask)

    #Create a boolean array based on indices of subbands channels
    arr_subband_index = subband_edges(bw, n_chan, tele)
    bool_array = np.zeros(ds_off.shape, dtype=bool)
    bool_array[arr_subband_index] = True

    #Mask channels where the value is zero and mask the channels where the boolean list is True
    ds_off_masked = np.ma.masked_where(ds_off == 0.00, ds_off)
    ds_off_masked = np.ma.masked_array(ds_off_masked, mask=bool_array)
    ds_corrected = (ds - ds_off_masked.mean(axis=1)[:,None]) / ds_off_masked.std(axis=1)[:, None]
    return ds_corrected


def SN(ts, t_cent, on_mask=None, t_sig=None):
    # duplicate the on-mask which comes as 2D but have it apply to a time series only
    on_mask_t = on_mask.T.mean(axis=-1)
    ts_off = np.ma.masked_array(ts, mask=on_mask_t)

    ts -= np.mean(ts_off)
    ts_off -= np.mean(ts_off)
    ts /= np.std(ts_off)
    ts_off /= np.std(ts_off)

    return ts, on_mask_t


def plot(ds, extent, tsamp, n_chan, tele, plot_grid_idx, fig, width,
         burst_n, nrows, ncols, nbursts, index,
         t_cent, f_cent,
         t_start, t_end,
         t_sig=None, f_sig=None,
         plot_spectrum=False, colour='cyan',
         t_cent_2=None,t_sig_2=None,f_cent_2=None,f_sig_2=None):
    """
    Creates the family burst plot
    """
    bw = extent[3]-extent[2]

    #conv = 2.355 #conversion from sigma to FWWHM
    #conv = 2 * np.sqrt(2 * (np.log(2)))


    on_mask = np.zeros_like(ds, dtype=bool)
    for s, e in zip(t_start, t_end):
        on_mask[:, s:e] = True

    #function to normalise/standardize the data
    ds_norm = normalise(ds, t_cent, bw, n_chan, tele, on_mask)

    # on_mask is True where the pulse in 'on', masked entries
    # are those where the mask is True. Thus the 'off' mask is
    # the inverse of on_maks
    ds_norm_on = np.ma.masked_array(ds_norm, mask=~on_mask)

    spectrum = np.mean(ds_norm_on, axis=-1)

    #Add the offset from SFXC
    res_f = (extent[3]-extent[2])/spectrum.size #frequency resolution

    if t_cent_2:
        t_cent_1 = np.min([t_cent,t_cent_2])
        if t_cent_2 == t_cent_1:
            t_cent_2 = t_cent
            t_sig_1 = t_sig_2
            t_sig_2 = t_sig
        else:
            t_sig_1 = t_sig

        t_cent = (t_cent_2-t_cent_1)/2. + t_cent_1

    if f_cent_2:
        f_cent_1 = np.min([f_cent,f_cent_2])
        if f_cent_2 == f_cent_1:
            f_cent_2 = f_cent
            f_sig_1 = f_sig_2
            f_sig_2 = f_sig
        else:
            f_sig_1 = f_sig
        f_cent = (f_cent_2-f_cent_1)/2.+ f_cent_1

    if f_cent and f_sig:
        f_l_bin = int(np.ceil((f_cent-2.*f_sig)))# - extent[2])/res_f))
        if f_sig_2:
            f_l_bin = int(np.ceil((f_cent_1-2.*f_sig_1)))

        if f_l_bin < 0:
            f_l_bin = 0
        f_h_bin = int(np.ceil((f_cent+2.*f_sig)))# - extent[2])/res_f))
        if f_sig_2:
            f_h_bin = int(np.ceil((f_cent_2+2.*f_sig_2)))
        if f_h_bin > spectrum.size:
            f_h_bin = spectrum.size
    else:
        f_l_bin = 0
        f_h_bin = spectrum.size

    #bw = (f_h_bin - f_l_bin)*res_f #MHz


    #Create time window for burst, redfine the extent
    peak = np.int(np.ceil(t_cent))
    extent[0] = - width / 2. # width is time window around the burst
    extent[1] = width / 2.

    #time bins of burst around window
    t_l_bin = int(t_cent - (width/2.)/(tsamp*1e3))
    t_h_bin = int(t_cent + (width/2.)/(tsamp*1e3))
    if t_l_bin < 0 or t_h_bin > ds_norm.shape[1]:
        print(f'\n\n Careful with the on-times, we are rolling.\n')
        ds_norm = np.roll(ds_norm, int(ds_norm.shape[1]/2.-t_cent), axis=1)
        on_mask = np.roll(on_mask, int(ds_norm.shape[1]/2.-t_cent), axis=1)

        t_start_offsets = t_start - t_cent
        t_end_offsets = t_end - t_cent

        t_cent = int(ds_norm.shape[1]/2.)

        t_start =  t_cent + t_start_offsets
        t_end = t_cent + t_end_offsets
        t_h_bin = int(t_cent + (width/2.)/(tsamp*1e3))
        t_l_bin = int(t_cent - (width/2.)/(tsamp*1e3))

    #Time series
    ts = np.mean(ds_norm[f_l_bin:f_h_bin,:], axis=0) # time series (summed the frequencies) around the burst
    ts, on_mask_t = SN(ts, t_cent, on_mask)
    ds_norm = ds_norm[:,t_l_bin:t_h_bin]
    on_mask = on_mask[:,t_l_bin:t_h_bin]
    ts = ts[t_l_bin:t_h_bin]
    print(f't_l_bin: {t_l_bin}; t_cent: {t_cent}; t_h_bin: {t_h_bin}')
    on_mask_t = on_mask_t[t_l_bin:t_h_bin]
    t_start -= int(t_cent)
    t_end -= int(t_cent)
    print(f'After fiddling: t_start = {t_start}\nt_end={t_end}')
    print(f'In time units: t_start = {t_start*(tsamp*1e3)}\nt_end={t_end*(tsamp*1e3)}')
    x = (np.arange(t_l_bin,t_h_bin,1)-int(t_cent))*tsamp*1e3
    print(f'x_min = {x.min()}\nx_max = {x.max()}')
    #Gridspec in Gridspec
    rows = 2
    cols = 1
    if plot_spectrum: cols += 1
    plot_grid = gridspec.GridSpecFromSubplotSpec(rows, cols, plot_grid_idx, wspace=0., hspace=0.,
                                                 height_ratios=[1,]*(rows-1)+[2,], width_ratios=[5,]+[1,]*(cols-1))
    ax1 = plt.Subplot(fig, plot_grid[rows-1,0])
    ax2 = plt.Subplot(fig, plot_grid[rows-2,0], sharex=ax1)
    if plot_spectrum: ax3 = plt.Subplot(fig, plot_grid[rows-1,1], sharey=ax1)
    units = ("GHz", "ms")

    ##Making the marker sizes
    #We need to copy the colormap, old command not usable anymore #cmap = plt.cm.viridis
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))

    #size of the red marker, 1 if full freq channel, 0 is no red mask -> just an int
    zap_size = int(ds_norm.shape[1]/18)

    #identifying the frequency channels that where zapped, only need a 1-D boolean list
    zapped = ma.getmaskarray(ds_norm[:,1])

    #We take the bins where the red marker will be and replace the mask with values of 9999
    #We do this because we want to asign a individual color map to these values
    zapped_marker = ds_norm[zapped,:zap_size]
    zapped_marker_mask_value = 1999
    ds_norm[zapped,:zap_size] = ma.filled(zapped_marker,fill_value=zapped_marker_mask_value)

    #Here we create the colormap of red
    cm1 = mpl.colors.ListedColormap(['red'])

    #This is a boolean 2d array wehre values of 9999 at True
    mask_marker = ds_norm==zapped_marker_mask_value

    #We now create a 2D masked array where the values are NOT 9999, so we have a masked array with only True's
    #This means that we only make the marker red, the white spaces in the plot at part of the original mask
    mask_marker = ma.masked_where(mask_marker==False,mask_marker)

    #setting up quantiles, using quantiles we can supress the super bright and dim parts of the dynamic spectrum
    #this way the constrast of the plot improves and the burst is better visible.
    #We first replace the masking with nans, because we can then use nanpercentile
    replace_masked = ds_norm.filled(np.nan)
    quantiles = np.nanpercentile(replace_masked, (2, 98))


    #exception for B07-o8 due to RFI
    if burst_n == "B07-o8":
        quantiles = np.nanpercentile(replace_masked, (2, 95))
        print("no quantiles used, this is low bit data")

    #Plotting the dynamic spectrum
    ax1.imshow(ds_norm, cmap=cmap, origin='lower', aspect='auto', interpolation='nearest',
               vmin=quantiles[0], vmax=quantiles[1],extent=extent)

    #plotting the red markers for the masked channels
    ax1.imshow(mask_marker, cmap=cm1, origin='lower', aspect='auto',                interpolation='nearest',                extent=extent)

    ax1.set_xlim(extent[0]-0.0001, extent[1]+0.0001)
    ax1.set_ylim(extent[2],extent[3])

    #Label only edge plots
    if index % ncols == 0:
        ax1.set_ylabel(r'${\rm Frequency}\ ({\rm GHz})$')
        ax2.set_ylabel('S/N')
        #ax1.set_yticklabels([r'$1.65$',r'$1.70$'])
    #else:
    #    ax1.tick_params(axis='y', labelleft=False)

    ax1.tick_params(axis='x', labelbottom=False)

    ax1.yaxis.set_major_locator(LinearLocator(numticks=4))
    y_ticks_ghz = ax1.get_yticks() / 1000.
    y_ticks_ghz_round = np.around(y_ticks_ghz, decimals=2)
    ax1.set_yticklabels(y_ticks_ghz_round)

    #Label only edge plots
    if index >= ncols * (nrows-1):
        ax1.tick_params(axis='x', labelbottom=True)
        ax1.set_xlabel(r'${\rm Time}\ ({\rm ms})$')

    #ax1.tick_params(axis='y',labeltop='off', labelleft='off')
    ax1.xaxis.set_minor_locator(MultipleLocator(5))

    ##Plotting SNR ax2
    #plot pulse profile (time series)
    ax2.plot(x, ts, 'k-',alpha=1.0,zorder=1,lw=0.5, drawstyle='steps-mid')
    point = ax2.scatter(x[0], ts[0], facecolors='none', edgecolors='none')

    #Setting ranges for SNR
    y_range = ts.max() - ts.min()
    ax2.set_ylim(-y_range/3., ts.max()*1.1)
    maxSN = np.ceil(np.max(ts))

    #only add 2 ticks
    yticks_sn=np.array([0,maxSN])

    ax2.set_yticks(yticks_sn)

    if t_sig_2:
        #ax2.hlines(y=-y_range/3.,xmin=(-(t_cent-t_cent_1+t_sig_1*2)*(tsamp*1e3)),xmax=((t_cent_2-t_cent+t_sig_2*2)*(tsamp*1e3)), lw=10,color=colour,zorder=0.8, alpha=0.2)
        ax2.hlines(y=-y_range/3.,xmin=(-(t_cent-t_cent_1+t_sig_1)*(tsamp*1e3)),xmax=((t_cent_2-t_cent+t_sig_2)*(tsamp*1e3)), lw=10,color=colour,zorder=0.8) #, alpha=0.4)
    else:
        #ax2.hlines(y=-y_range/3.,xmin=(-t_sig*2*(tsamp*1e3)),xmax=(t_sig*2*(tsamp*1e3)), lw=10,color=colour,zorder=0.8, alpha=0.2)
        ax2.hlines(y=[-y_range/3.]*len(t_start),
                   xmin=(t_start*(tsamp*1e3)), xmax=(t_end*(tsamp*1e3)),
                   lw=10, color=colour, zorder=0.8) #, alpha=0.4)

    #b = np.argmin(np.abs(x-(-t_sig*2*(tsamp*1e3))))
    #e = np.argmin(np.abs(x-(t_sig*2*(tsamp*1e3))))
    if t_sig_2:
        b = np.argmin(np.abs(x-(-((t_cent-t_cent_1)+t_sig_1*2)*(tsamp*1e3))))
        e = np.argmin(np.abs(x-((t_cent_2-t_cent+t_sig_2*2)*(tsamp*1e3))))

    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width*0.65, box.height])
    legend1=ax2.legend((point,point), ((r"${\rm %s}$")%burst_n,""),
                       loc='upper left',handlelength=0,bbox_to_anchor=(0.01, 1.05),
                       handletextpad=-0.5,frameon=False,markerscale=0, fontsize = 6)
    ax2.add_artist(legend1)
    ax2.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                    bottom=True, top=False, left=True, right=False)

    legend2=ax2.legend((point,point), ((fr"{int(tsamp*1e6)} $\mu$s"),""),
                       loc='upper left',handlelength=0,bbox_to_anchor=(0.6, 1.05),
                       handletextpad=-0.5,frameon=False,markerscale=0,fontsize=6)
    ax2.add_artist(legend2)

    legend3=ax2.legend((point,point), ((f"{res_f:.2f} MHz"),""), loc='upper left',
                       handlelength=0,bbox_to_anchor=(0.6, 0.9), handletextpad=-0.5,
                       frameon=False,markerscale=0,fontsize=6)
    ax2.add_artist(legend3)

    ax2.tick_params(axis='y', which='major', pad=1.5)

    ##Plotting spectrum ax3
    #plot spectrum (amplitude vs freq) only if plot_spectrum=True
    if plot_spectrum:
        y = np.linspace(extent[2], extent[3], spectrum.size)
        ax3.plot(spectrum, y, 'k-',zorder=2, lw=0.7,drawstyle='steps-mid')
        ax3.set_ylim(extent[2],extent[3])
        x_range = spectrum.max() - spectrum.min()
        ax3.set_xlim(-x_range/3., x_range*6./5.)

        #Setting up the ticks
        ax3.tick_params(labelleft=False, labelbottom=False, labeltop=False, labelright=False,                        bottom=False, top=False, left=True, right=False)
        ax3.tick_params(which="both", direction='in')
    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    if plot_spectrum: fig.add_subplot(ax3)

    return


def prep_data(dataframe, mode='sfxc'):
    known_modes = ['sfxc', 'scale', 'spc']
    if not mode in known_modes:
        raise ValueError(f'Unknown mode {mode}. Chose from {known_modes}.')

    data = pd.read_csv(dataframe)
    data = data[(data.src == mode)]
    data.sort_values(by='toa_bary_tdb_inf_freq', ignore_index=True, inplace=True)
    ids = data.id.unique()
    return ids, data


if __name__ == '__main__':
    #adapt the figure parameters as needed:
    width = 80. #Time window around the burst in ms.

    nrows = 5
    ncols = 4

    #grid of burst plots
    #wspace and hspace is the space between plots
    plot_grid = gridspec.GridSpec(nrows, ncols, wspace=0.3, hspace=0.1)

    # base directory
    base_dir = '/data1/franz/'

    SEFD = 24./1.54
    dist = 453.

    ids, data = prep_data('../dbs/burst_info.csv')
    idx = 0
    part = 0
    nplots = 0
    #colours = ['darkorchid', '#377eb8', '#4dac26', '#ca0020'] # [st, wb, o8, tr]
    clrs = get_station_colors()
    colours = [clrs['st'], clrs['wb'], clrs['on'], clrs['tr']]
    dishes = ['st', 'wb', 'o8', 'tr']

    # for better visibility need to average down in freq and/or time for some bursts
    tscrunch_ids = {'B02-o8': 2,
                    'B03-o8': 2,
                    'B04-o8': 4,
                    'B07-o8': 8,
                    'B09-o8': 4,
                    'B11-o8': 4,
                    'B14-o8': 8,
                    'B15-wb': 8,
                    'B16-wb': 8,
                    'B17-wb': 2,
                    'B21-wb': 4,
                    'B23-wb': 2,
                    'B25-wb': 2,
                    'B26-o8': 2,
                    'B27-o8': 2,
                    'B26-wb': 2,
                    'B28-wb': 2,
                    'B30-wb': 2,
                    'B32-o8': 2,
                    'B33-st': 2,
                    'B33-wb': 2,
                    'B35-wb': 2,
                    'B35-st': 2,
                    'B36-wb': 2,
                    'B38-wb': 4,
                    'B39-st': 4,
                    'B39-wb': 2,
                    'B40-st': 2,
                    'B41-wb': 2,
                    'B46-st': 2
                    }

    fscrunch_ids = {'B02-o8': 8,
                    'B03-o8': 8,
                    'B04-o8': 8,
                    'B07-o8': 16,
                    'B09-o8': 8,
                    'B11-o8': 8,
                    'B14-o8': 8,
                    'B15-wb': 8,
                    'B16-wb': 8,
                    'B17-wb': 4,
                    'B21-wb': 16,
                    'B23-wb': 8,
                    'B24-wb': 8,
                    'B25-wb': 16,
                    'B26-o8': 16,
                    'B26-wb': 16,
                    'B27-o8': 8,
                    'B27-wb': 4,
                    'B28-wb': 16,
                    'B30-wb': 16,
                    'B32-o8': 8,
                    'B33-wb': 8,
                    'B35-wb': 4,
                    'B36-wb': 16,
                    'B39-wb': 4,
                    'B41-wb': 4,
                    'B44-st': 2
                    }

    # not sure why but some of the archives have the dedispersion already applied
    dm_zero = ['B02-o8', 'B03-o8', 'B07-o8', 'B14-o8']

    #size of the plot; for paper format (enough room for comments)
    fig = plt.figure(figsize=[10,12])

    for burst in tqdm(ids):
        #if not burst == 'B07-o8':
        #    continue
        print("--------------------")
        print(f'{burst}')
        print("--------------------")
        nplots += 1
        df = data[(data.id == burst)]
        experiment = df.experiment.values[0]
        dish = df.dish.values[0]
        scan = df.scan.values[0]
        colour = colours[dishes.index(dish)]
        tscrunch = None
        fscrunch = 2
        dm = 410.8
        if burst in tscrunch_ids.keys():
            tscrunch = tscrunch_ids[burst]
        if burst in fscrunch_ids.keys():
            fscrunch = fscrunch_ids[burst]
        if burst in dm_zero:
            dm = 0.0
        ds, extent, tsamp, n_chan = load_archive(f'{base_dir}/{experiment}/sfxc/{experiment}_{dish}_{scan}',
                                                 dm=dm, remove_baseline=False,
                                                 extent=True,
                                                 tscrunch=tscrunch,
                                                 fscrunch=fscrunch)
        print(tsamp, ds.shape)
        #create StokesI
        #shape is 3 for full pol and 2 for stokesI
        if len(ds.shape) > 2:
            ds_I = ds[0,:,:]+ds[1,:,:]

            print("extent before", extent)
            #Correct for SFXC bug
            extent[2] = extent[2] - 0.0625
            extent[3] = extent[3] - 0.0625
            print("extent after, ", extent, " ,we do this to account for an offset that SFXC introduces")

        else:
            print('Stockert burst')
            ds_I = ds

        # determine location of the burst in the file,
        # find the middle between components and convert to bin-space
        t0 = df.t0_ms.values
        print(f't0 = {t0} ms')
        t0 = (t0.min() + t0.max()) / 2. / (tsamp*1e3)
        print(f't0 = {t0} bins')
        # figure out the 'on'-time for each component from begin/end_acf_range
        # that one also comes in ms into the file, i.e. convert to bin-space too
        t_start = (df.begin_acf_range.values / (tsamp*1e3)).astype(int)
        t_end = (df.end_acf_range.values / (tsamp*1e3)).astype(int)
        print(f't_start = {t_start}; \nt_end = {t_end}')

        # we plot the specrum as the sum of all components
        # we don't indicate an 'on'-frequency range
        f0 = n_chan / 2

        plot(ds=ds_I, extent=extent, tsamp=tsamp, n_chan=n_chan,
             tele=dish, plot_grid_idx=plot_grid[idx], fig=fig, width=width,
             burst_n=burst, nrows=nrows, ncols=ncols, nbursts=len(ids), index=idx,
             t_cent=t0, f_cent=f0,
             t_start=t_start, t_end=t_end,
             plot_spectrum=True,
             colour=colour)

        idx+=1
        if (idx == nrows * ncols) or (nplots == len(ids)):
            idx = 0
            #left, right, bottom and top are used for space around the plot.
            fig.subplots_adjust(hspace=0.8, wspace=0.5, left=0.07,right=.98,bottom=.055,top=.99)

            fig.savefig(f"../plots/burst_family_pt{part}.pdf", format = 'pdf', dpi=150)
            fig.savefig(f"../plots/burst_family_pt{part}.png", facecolor='white', transparent=False,format = 'png')
            try:
                fig.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/burst_family_pt{part}.pdf",
                             format = 'pdf', dpi=150)
                fig.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/burst_family_pt{part}.png",
                            facecolor='white', transparent=False,format = 'png')
            except:
                print('Could not write to overleaf repo. Saved a local copy only.')
            part += 1

            #size of the plot; for paper format (enough room for comments)
            fig = plt.figure(figsize=[10,12])
