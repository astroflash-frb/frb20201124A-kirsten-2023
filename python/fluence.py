from radiometer import radiometer
from load_file import load_archive, load_filterbank
from interactive_sub_burst_identifier_1d import identify_bursts
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


def options():
    parser = argparse.ArgumentParser()
    general = parser.add_argument_group('General info about the data.')
    general.add_argument('arfile', type=str, default=None,
                         help='Archive file to perform the ACF analysis on')
    general.add_argument('-d', '--dm', type=float, default=0,
                         help='DM to dedisperse the data to.')
    general.add_argument('--tscrunch', type=int, default=1,
                         help='Downsampling factor in time for archive files.' + \
                         'Default=%(default)s')
    general.add_argument('--fscrunch', type=int, default=1,
                         help='Downsampling factor in frequency for archive files.' + \
                         'Default=%(default)s')
    general.add_argument('-t', '--begintime', type=int, default=0,
                         help='begin time bin to begin chopping in time')
    general.add_argument('-T', '--endtime', type=int, default=0,
                         help='end time bin to end chopping in time')
    general.add_argument('-b', '--beginchan', type=int, default=0,
                         help='begin freq bin to begin chopping in freq')
    general.add_argument('-e', '--endchan', type=int, default=0,
                         help='end freq bin to end chopping in freq')
    general.add_argument('-D', '--distance', type=float, default=None,
                         help='Distance to the source in Mpc')
    general.add_argument('--filterbank',action='store_true',
                         help='if input is in filterbank format, give this option [default: False, assuming archive file input].')
    general.add_argument('-S', '--sefd', type=float, default=20/1.54,
                         help='system equivalent flux density [default: 20/1.54 for Eff]')
    parser.add_argument('--db', type=str, default=None,
                        help='Path to pickle file that contains the pandas data frame.')
    parser.add_argument('--src', type=str, default='any',
                        help='Used when writing to a pandas data frame. In case parameters ' +
                        'for the same burst are determined in different ways (e.g. different files)' +
                        'this parameter can be used to distinguish them.')
    return parser.parse_args()


def fluences(ts, tsamp, bw, SEFD, distance=None):
    """
    ts is already calibrated time series -- i.e. in S/N units
    already cropped in time and frequency

    """
    profile_flux = ts*radiometer(tsamp*1000, bw, 2, SEFD)
    fluence = np.sum(profile_flux*tsamp*1000)
    peak_snr = np.max(ts)
    peak_flux = np.max(profile_flux)

    print("Peak S/N is "+str(np.max(ts)))
    print("Peak flux is "+str(np.max(profile_flux))+" Jy")
    print("Fluence is "+str(fluence)+" Jy ms")

    if distance is not None:
        #convert Jy ms to J s
        fluence_Jys = fluence*1e-3
        #convert Mpc to cm
        distance_lum_cm = 3.086e24*distance
        energy_iso= fluence_Jys*4*np.pi*(distance_lum_cm**2)*1e-23
        lum_spec = energy_iso/(len(profile_flux)*tsamp)
        print("Isotropic energy is "+str(energy_iso)+" erg Hz^-1")
        print("Spectral luminosity is "+str(lum_spec)+" erg s^-1 Hz^-1")
        return fluence, peak_snr, peak_flux, energy_iso, lum_spec

    return fluence, peak_snr, peak_flux


if __name__ == "__main__":
    args = options()

    # load in archive
    if args.filterbank is False:
        archive, extent, tsamp, T0 = load_archive(args.arfile, dm=args.dm,
                                                  extent=True,
                                                  fscrunch=args.fscrunch,
                                                  tscrunch=args.tscrunch)
        t_min, t_max, f_min, f_max = extent
        if f_min > f_max:
            f = f_min
            f_min = f_max
            f_max = f
        fsamp = (f_max - f_min) / archive.shape[0]  # in MHz
        # tsamp comes in seconds
    else:
        archive, extent, tsamp = load_filterbank(args.arfile, dm=args.dm,
                                                 fullpol=False)

    if len(archive.shape) == 3:
        archive = archive[0, :, :]+archive[1, :, :]  # only interested in Stokes I
    if args.db is not None:
        src = args.src
        data = pd.read_pickle(args.db)
        try:
            exp, dish, scan = (os.path.basename(args.arfile)).split('_')[0:3]
        except:
            raise ValueError('File name does not agree with assumed convention.')
        df = data[(data.experiment == exp) &
                  (data.dish == dish) &
                  (data.scan == scan) &
                  (data.src == src)]
        if not df.empty:
            central_t_bins = np.array(df.t0_ms.values / (tsamp * 1e3))
            central_f_bins = (df.f0_MHz.values - f_min) / fsamp
            off_range_t = np.array([df.off_range_t0.values[0],
                                    df.off_range_t1.values[0]]) // args.tscrunch

            hpbw_t_bins = (2.355 * np.abs(df.twidth_sigma_ms.values)) // (tsamp*1e3) + 1
            hpbw_f_bins = (2.355 * np.abs(df.fwidth_sigma_MHz.values)) // fsamp + 1
            # for the width in time and frequency we use the twice the HPBW
            ranges_t = np.array([[int(t_bin - hpbw_t),
                                  (t_bin + hpbw_t+1)] for t_bin, hpbw_t in zip(central_t_bins, hpbw_t_bins)],
                                dtype=np.int)
            # for multi-component burst, to make sure there is no
            # overlap between the ranges, i.e. to not count
            # fluences twice, we set the edges to not overlap
            if len(ranges_t) > 1:
                for i in range(len(ranges_t)):
                    end = ranges_t[i][1]
                    begin = ranges_t[i+1][0]
                    if begin < end:
                        middle = () // 2
                        ranges_t[i][1] = middle
                        ranges_t[i+1][0] = middle
                    if i+1 == len(ranges_t)-1:
                        break
            ranges_f = np.array([[f_bin - hpbw_f,
                                  f_bin + hpbw_f] for f_bin, hpbw_f in zip(central_f_bins, hpbw_f_bins)],
                                dtype=np.int)
            if off_range_t[0] == off_range_t[1] == 0:
                print(f"Please choose the off-pulse region (in time).")
                subbursts = identify_bursts(np.mean(archive, axis=0),
                                            "Please choose the off-pulse region (in time).")
                off_range_t = np.array(subbursts.ranges, dtype=np.int)
        else:
            print(f'Found no entry for {exp}_{dish}_{scan}. Doing it the manual way.')
    else:
        begintime = args.begintime
        if args.endtime != 0:
            if args.endtime <= args.begintime:
                raise ValueError('Please give an end time after the begin time, i.e. -T > -t')
            else:
                endtime = args.endtime
            ranges_t = np.array([begintime, endtime], dtype=np.int).reshape(-1, 2)
        else:
            subbursts = identify_bursts(np.mean(archive, axis=0))
            ranges_t = np.array(subbursts.ranges, dtype=np.int).reshape(-1, 2)

        print(f"Please choose the off-pulse region (in time).")
        subbursts = identify_bursts(np.mean(archive, axis=0))
        off_range_t = np.array(subbursts.ranges, dtype=np.int)

        beginchan = args.beginchan
        if args.endchan != 0:
            if args.endchan <= args.beginchan:
                raise ValueError('Please give an end frequency after the begin frequency, i.e. -e > -b')
            else:
                endchan = args.endchan
            ranges_f = np.array([beginchan, endchan], dtype=np.int).reshape(-1, 2)
        else:
            ranges_f = np.zeros_like(ranges_t)
            for component, Range in enumerate(ranges_t):
                subbursts = identify_bursts(np.mean(archive[:, Range[0]:Range[1]], axis=1))
                ranges_f[component] = np.array(subbursts.ranges, dtype=np.int)
    archive_origin = archive.copy()
    for component, Range in enumerate(ranges_t):
        archive = archive_origin.copy()
        beginchan, endchan = ranges_f[component]
        begintime, endtime = Range
        if beginchan < 0:
            beginchan = 0
        if begintime < 0:
            begintime = 0
        if endchan > archive.shape[0]:
            endchan = archive.shape[0] - 1
        if endtime > archive.shape[1]:
            endtime = archive.shape[1] - 1

        print(beginchan, endchan, begintime, endtime)
        archive_orig = archive[beginchan:endchan, :].copy()
        archive = archive[beginchan:endchan, begintime:endtime]

        profile = np.mean(archive, axis=0)
        profile_full = np.mean(archive_orig, axis=0)

        profile_off = np.mean(archive_orig[:, off_range_t[0]:off_range_t[1]],
                              axis=0)

        #convert to S/N
        profile -= np.mean(profile_off)
        profile_full -= np.mean(profile_off)
        profile_off -= np.mean(profile_off)
        profile /= np.std(profile_off)
        profile_full /= np.std(profile_off)

        plt.plot(profile_full)
        plt.ylabel('S/N')
        plt.axvline(off_range_t[0], color='red', label='off-region')
        plt.axvline(off_range_t[1], color='red')
        plt.axvline(begintime, color='blue', label='on-region')
        plt.axvline(endtime, color='blue')
        plt.legend()
        plt.show()
        plt.clf()
        vmin = 0
        vmax = archive_origin.var() * 6
        plt.imshow(archive_origin, aspect='auto', interpolation='Gaussian',
                   origin='lower', extent=extent, vmin=vmin, vmax=vmax)
        plt.axhline(f_min + beginchan*fsamp, color='green')
        plt.axhline(f_min + endchan*fsamp, color='green')
        plt.axvline((t_min + off_range_t[0]*tsamp) * 1000, color='red', label='off-region')
        plt.axvline((t_min + off_range_t[1]*tsamp) * 1000, color='red')
        plt.axvline((t_min + begintime*tsamp) * 1000, color='red', label='on-region')
        plt.axvline((t_min + endtime*tsamp) * 1000, color='red')
        plt.show()
        plt.clf()
        #convert to Jy units
        nchans = archive.shape[0]
        bw = nchans * fsamp
        SEFD = args.sefd

        try:
            vals = fluences(profile, tsamp, bw, SEFD, distance=args.distance)
            fluence, peak_snr, peak_flux, energy_iso, lum_spec = vals
        except:
            vals = fluences(profile, tsamp, bw, SEFD, distance=args.distance)
            fluence, peak_snr, peak_flux = vals

        if args.db is not None:
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'fluence_jyms'] = fluence
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'peak_snr'] = peak_snr
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'peak_flux_jy'] = peak_flux
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'energy_iso_erg_per_hz'] = energy_iso
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'spectral_lum_erg_per_s_per_hz'] = lum_spec
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'off_range_t0'] = off_range_t[0] * args.tscrunch
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.src == src) &
                     (data.scan == scan) & (data.component == component),
                     'off_range_t1'] = off_range_t[1] * args.tscrunch
    if args.db is not None:
        data.to_pickle(args.db)
