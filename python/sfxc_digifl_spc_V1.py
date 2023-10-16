#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spopt
import numpy.ma as ma
from tqdm import tqdm
from load_file_kenzie import load_archive
from matplotlib import transforms
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
import os

rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['axes.labelsize'] = 12
rcParams['legend.fontsize'] = 10
rcParams['legend.handlelength'] = 1.0


def subband_edges(bandwidth, n_chan):
    #[   0    1  126  127  128  129  254  255  256  257  382  383  384  385
    #510  511  512  513  638  639  640  641  766  767  768  769  894  895
    #896  897 1022 1023]
    # ab is starting array [0,1,126,127]
    subband_edge = n_chan // 8
    #if statement so to make sure then bandpass also gets smaller when downsampling
    ab = np.array([0,1,subband_edge-2,subband_edge-1])

    if n_chan == 512:
        ab = np.array([0,subband_edge-1])
    if n_chan <= 256:
        ab = np.array([0])
    print("bandwidth is:", bandwidth, ". Number of channels:", n_chan, ". First subband array: ", ab)
    print("Subband edges: ", subband_edge)

    mask = []
    for i in range(8):
        a = ab+(i*subband_edge)
        mask.append(a)
    mask_arr = np.asarray(mask).ravel()
    return mask_arr


def correcting_data(ds_burst, n_chan,
                    downsample, bandwidth,
                    region_1, region_2, mask=None):
    #bins where there is not burst present
    ds_woburst_1 = ds_burst[:, region_1[0]//downsample:region_1[1]//downsample]
    ds_woburst_2 = ds_burst[:, region_2[0]//downsample:region_2[1]//downsample]
    ds_woburst = np.concatenate((ds_woburst_1, ds_woburst_2), axis=1)

    #make indices to mask the subbands
    arr_subband_index = subband_edges(bandwidth, n_chan)

    #Create a boolean array based on indices of subbands channels
    bool_array = np.zeros(ds_woburst.shape, dtype=bool)
    bool_array[arr_subband_index] = True
    if mask:
        #load the mask file in, take the first element so it is only the array
        #flip it, because shit is flipped
        b_array_mask = np.flip(np.load(mask)[0])
        bool_array[b_array_mask == 0] = True
    #Mask channels where the value is zero and mask the channels where the boolean list is True
    ds_woburst_masked = ma.masked_where(ds_woburst == 0.00, ds_woburst)
    ds_woburst_masked = ma.array(ds_woburst_masked, mask=bool_array)

    #Calculate the mean and the standard deviation
    #apply the correction -> with mean and standard deviation
    ds_corrected = (ds_burst - ds_woburst_masked.mean(axis=1)[:,None]) / ds_woburst_masked.std(axis=1)[:, None]
    return ds_corrected


def sum_loop(file, downsample, f_fscrunch, region_1, region_2, bandwidth,
             maskfile=None, dm_value=None):
    #load in the archive file with a certain DM value
    if dm_value == None:
        arr, extent, tsamp, n_chan = load_archive(file, tscrunch=downsample, fscrunch=f_fscrunch,
                                                  extent=True)
    else:
        arr, extent, tsamp, n_chan = load_archive(file, tscrunch=downsample, fscrunch=f_fscrunch,
                                                  dm=dm_value, extent=True)
    if len(arr.shape) < 3:
        print("Stokes I archive file")
        #StokesI archive file
        I = arr
    else:
        print("extent before", extent)
        #Correct for SFXC bug
        extent[2] = extent[2] - 0.0625
        extent[3] = extent[3] - 0.0625
        print("extent after, ", extent, " ,we do this to account for an offset that SFXC introduces")
        #create StokesI
        I = arr[0,:,:]+arr[1,:,:]

    #Downsampling
    #I_down = I.reshape(n_chan,-1,downsample).mean(axis=-1)
    #Correct for the baseline, calculate statistics
    if maskfile == None:
        I_br_down = correcting_data(I, n_chan, downsample, bandwidth, region_1, region_2, mask=None)
    else:
        I_br_down = correcting_data(I, n_chan, downsample, bandwidth, region_1, region_2, mask=maskfile)
    return I_br_down, extent, tsamp, n_chan


def plot_corrected2(ds_corrected, extent, tsamp, n_chan,
                    downsample, fscrunch, region_1, region_2, bandwidth,
                    plt_save,
                    xlim_l=None, xlim_r=None, ylim_b=None, ylim_t=None, roll=None):
    #Calculate the quantiles, take into account the masking of the array
    replace_masked = ds_corrected.filled(np.nan)
    new_quantiles = np.nanpercentile(replace_masked, (1, 99))

    #information on time resolution, and downsampling factor in us
    #tsamp_d is in miliseconds
    tsamp_d = tsamp * 1000.0
    #tsamp_d = downsample * tsamp * 1000.0
    #width in KHz
    channel_width = (bandwidth/n_chan)
    #roll = 15000
    freqs = np.linspace(extent[2],extent[3], endpoint=False, num=ds_corrected.shape[0])
    #Freq_channels, based on extent
    #sum_freq = ma.MaskedArray.sum(ds_corrected, axis=0)

    if ylim_b == None:
        sum_freq = ma.MaskedArray.sum(ds_corrected, axis=0)
    else:
        #Find indices in y_range based no the limits set by user, then only take those ranges
        #sum_freq is now limited in frequency
        freq_index = np.where( (ylim_b <= freqs) & (freqs <= ylim_t) )
        sum_freq = ma.MaskedArray.sum(ds_corrected[freq_index], axis=0)

    #Take the regions without burst and add them together
    sum_freq_woburst_1 = sum_freq[region_1[0]//downsample:region_1[1]//downsample]
    sum_freq_woburst_2 = sum_freq[region_2[0]//downsample:region_2[1]//downsample]
    sum_freq_woburst = np.concatenate((sum_freq_woburst_1, sum_freq_woburst_2))

    #calculate the snr of the frequency range
    snr = (sum_freq-sum_freq_woburst.mean())/sum_freq_woburst.std()

    #Create x grid to plot the step function.
    x_grid = np.linspace(extent[0], extent[1], endpoint=False, num=ds_corrected.shape[1])

    if xlim_l == None:
        sum_time = ds_corrected.sum(axis=1)
    else:
        #Get the time bins where you want to zoom in
        time_index = np.where( (xlim_l <= x_grid) & (x_grid <= xlim_r))

        #recreate what the user sees, so you select the correct bins
        plotted_ds = np.roll(ds_corrected,roll//downsample, axis=1)

        #Now select the correct time bins
        sum_time_index = plotted_ds[:, time_index[0]]

        #new values for sum_time
        sum_time = sum_time_index.sum(axis=1)

    #Take the regions without burst and add them together
    time_woburst_1 = ds_corrected[: , region_1[0]//downsample:region_1[1]//downsample]
    time_woburst_2 = ds_corrected[: , region_2[0]//downsample:region_2[1]//downsample]

    tsamp_plot = tsamp * 1e6
    #Axis=1 because of 2D array
    time_woburst = np.concatenate((time_woburst_1, time_woburst_2), axis=1)
    sum_time_snr = (sum_time - time_woburst.mean(axis=1)) / time_woburst.std(axis=1)
    fig = plt.figure(figsize=(5, 5))
    #widths = [15, 2]
    heights = [1, 4]
    spec5 = fig.add_gridspec(ncols=1, nrows=2,
                             height_ratios=heights,
                             wspace=0.0, hspace=0.0)
    ax1 = fig.add_subplot(spec5[0, 0])
    #ax1.set_title(f"Dynamic spectrum of R67 burst pre028-wb \n DM = 410.775 -- {tsamp_plot:.2f} $\mu$s -- {channel_width} MHz channels",                  fontsize=22, pad=20)
    ax1.set_ylabel('S/N')
    ax1.tick_params()
    ax1.plot(x_grid, np.roll(snr, roll//downsample),
             drawstyle='steps-pre',color="black")
    ax1.tick_params(labelbottom=False)
    ax1.set_ylim(-10, 110)
    ax2 = fig.add_subplot(spec5[1, 0], sharex=ax1)
    ax2.imshow(np.roll(ds_corrected, roll//downsample, axis=1),
               aspect="auto", interpolation="none",
               vmin=new_quantiles[0], vmax=new_quantiles[1],origin='lower',
               extent=(extent))
    ax2.set_ylabel('Frequency (MHz)')
    ax2.set_xlabel('Time (ms)')
    #ax2.tick_params(labelsize=14)
    #ax2.tick_params(labelsize=14)
    ax2.set_xlim(xlim_l, xlim_r)
    ax2.set_ylim(ylim_b, ylim_t)
    middle = (xlim_l + xlim_r) // 2
    nticks = 7
    xtick_lables = np.linspace(-30, 30, nticks, dtype=int).astype(str)
    xticks = np.linspace(middle-30, middle+30, nticks)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_lables)

    #ax3 = fig.add_subplot(spec5[1, 1],sharey=ax2)
    #ax3.plot(sum_time_snr, freqs, drawstyle='steps-pre', color="black")
    #ax3.tick_params(labelbottom=False, labelleft=False)

    plot_name = plt_save + f"_{tsamp_plot:.2f}us" + f"_{channel_width}MHz"
    plotdir = '../plots/'
    plt.savefig(plotdir + plot_name + ".pdf",  bbox_inches='tight',
                facecolor='white', transparent=False, dpi=600)
    try:
        plotdir = f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/"
        plt.savefig(plotdir + plot_name + ".pdf",  bbox_inches='tight',
                    facecolor='white', transparent=False, dpi=600)
    except:
        print('Could not save plot to overleaf repo. Got a local copy only.')
    #plt.savefig(plot_name + ".jpeg")
    return snr, time_woburst, fig


def archive_dif(snr_scale, snr_scale_spc, snr_sfxc, show=False):
    max_index_scale = np.argmax(snr_scale)
    max_index_sfxc = np.argmax(snr_sfxc)
    index_diff = max_index_scale - max_index_sfxc
    zeros = np.zeros(index_diff)
    snr_sfxc_shifted = np.concatenate((zeros, snr_sfxc))
    len_diff = len(snr_scale) - len(snr_sfxc_shifted)
    snr_sfxc_shifted = np.concatenate((snr_sfxc_shifted, np.zeros(len_diff)))
    if show:
        plt.figure(figsize=(16, 8))
        plt.plot(snr_scale, label="scale")
        plt.plot(snr_scale_spc, label="scale_spc")
        plt.plot(snr_sfxc_shifted,label="snr_sfxc")
        print(np.argmax(snr_sfxc))
        plt.legend()
        plt.xlim(7400,8100)
        plt.axhline(0)
        plt.ylabel("SNR")
        plt.xlabel("bins")
        plt.show()

    return snr_sfxc_shifted


def archive_dif_2(snr_scale, snr_scale_spc, snr_sfxc_shifted):
    #shift the x-axis
    x_range = np.arange(0, len(snr_scale))*128e-3
    x_range_shift = x_range - (0.991*1000)

    fig = plt.figure(figsize=(5, 5)) #, constrained_layout=True)
    spec = fig.add_gridspec(ncols=1, nrows=1)

    ax0 = fig.add_subplot(spec[0, 0])
    ax0.plot(x_range_shift, snr_scale, label="Digifil", color="#66c2a5", linewidth=1.5)
    ax0.plot(x_range_shift, snr_scale_spc, label="SPC", color="#fc8d62", linewidth=1.5)
    ax0.plot(x_range_shift, snr_sfxc_shifted, label="SFXC", color="#8da0cb", linewidth=1.5)
    ax0.set_ylabel("Signal-to-noise")
    ax0.axhline(0, color='black', linestyle='dashdot')
    ax0.set_xlim(-20,20)
    ax0.set_ylim(-8,105)
    ax0.legend()
    #ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Time [ms]")
    ax0.tick_params('both', length=5, width=1, which='major')
    plt.savefig("../plots/sfxc_digi_spc.png", bbox_inches='tight',
                facecolor='white', transparent=False, dpi=600)
    try:
        plt.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/sfxc_digi_spc.png",
                    bbox_inches='tight',
                    facecolor='white', transparent=False, dpi=600)
    except:
        print('Could not save plot to overleaf repo. Got a local copy only.')

    return fig




#scaled Digifil data, same params for both scale and spc
downsample = 1
fscrunch = 1
bandwidth = 128
region_1 = [1000, 5000]
region_2 = [12000, 15000]

I_br_down_sc, extent_sc, tsamp_sc, n_chan_sc = \
    sum_loop("../data/r67028_o8_no0020_allIFs.vdif_pol2.fil.scale.ds.ar",
             downsample, fscrunch, region_1, region_2, bandwidth,
             dm_value=410.775)

snr_scale, time_woburst_sc, fig_sc = \
    plot_corrected2(I_br_down_sc, extent_sc, tsamp_sc, n_chan_sc,
                    downsample, fscrunch, region_1, region_2, bandwidth,
                    plt_save="r67028_scale",
                    roll=0, xlim_l=956, xlim_r=1026, ylim_b=None, ylim_t=None)

I_br_down_spc, extent_spc, tsamp_spc, n_chan_spc = \
    sum_loop("../data/r67028_o8_no0020_allIFs.vdif_pol2.fil.scale.spc.ds.ar",
             downsample, fscrunch, region_1, region_2, bandwidth,
             dm_value=410.775)

snr_scale_spc, time_woburst_spc, fig_spc = \
    plot_corrected2(I_br_down_spc, extent_spc, tsamp_spc, n_chan_spc,
                    downsample, fscrunch, region_1, region_2, bandwidth,
                    plt_save="r67028_spc",
                    roll=0, xlim_l=956, xlim_r=1026, ylim_b=None, ylim_t=None)

#scaled SFXC data, different regions plotted,
# factor 16 accounts for the fact the regions were chosen
# on a file with 16 times the current time resolution
region_1 = [1000//16, 5000//16]
region_2 = [22000//16, 30000//16]

I_br_down_sf, extent_sf, tsamp_sf, n_chan_sf = \
    sum_loop("../data/r67028_o8_no0020_allIFs.sfxc.ds.ar",
             downsample, fscrunch, region_1, region_2, bandwidth,
             dm_value=410.775)

snr_sfxc, time_woburst_sfxc, fig_sfxc = \
    plot_corrected2(I_br_down_sf, extent_sf, tsamp_sf, n_chan_sf,
                    downsample, fscrunch, region_1, region_2, bandwidth,
                    plt_save="r67028_sfxc",
                    roll=0, xlim_l=61, xlim_r=131, ylim_b=None, ylim_t=None)

snr_sfxc_shifted = archive_dif(snr_scale, snr_scale_spc, snr_sfxc)
fig_comparison = archive_dif_2(snr_scale, snr_scale_spc, snr_sfxc_shifted)
