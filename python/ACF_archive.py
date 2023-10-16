
"""
Takes an archive file as input.
Computes the 2D ACF of the dynamic spectrum with 1D ACF plotted either side.
Tries to fit a 2D gaussian to the 2D ACF, removes this from the 1D ACFs and fits a Lorentzian to the residuals

Since it involves fitting 2D Gaussian to the 2D ACF and Lorentzians to the 1D residuals, you may need to change some of the initial guesses and ranges to use for the fits, depending on your situation.

Kenzie Nimmo 2021
"""

from ACF_funcs import autocorr_fft, autocorr_2D, lorentz
from load_file import load_archive
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import pandas as pd
from lmfit import minimize, Parameters, fit_report, Model
import os
from interactive_sub_burst_identifier_1d import identify_bursts
from interactive_sub_burst_identifier_2d import identify_bursts as identify_bursts_2d
from interactive_sub_burst_identifier_2d import ds
from init_db import init_pandas_df
from radiometer import radiometer


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


def gaus(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def twoD_Gaussian(x_data_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta):
    (x,y) = x_data_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()


def twoD_Gaussian_fit(params,x_data_tuple,data):
    amplitude=params['amplitude']
    xo=params['xo']
    yo=params['yo']
    sigma_x = params['sigma_x']
    sigma_y = params['sigma_y']
    theta = params['theta']

    fit=twoD_Gaussian(x_data_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta)

    resid = data.ravel()-fit
    return resid


def options():
    parser = argparse.ArgumentParser()
    general = parser.add_argument_group('General info about the data.')
    general.add_argument('arfile', type=str, default=None,
                         help='Archive file to perform the ACF analysis on')
    general.add_argument('-d','--dm',type=float, default=0,
                         help='DM to dedisperse the data to.')
    general.add_argument('-t', '--begintime', type=int, default=0,
                         help='begin time bin to begin chopping in time')
    general.add_argument('-T', '--endtime', type=int, default=0,
                         help='end time bin to end chopping in time')
    general.add_argument('-l', '--loadnpy',type=str,default=None,
                         help='Load in a numpy file containing the ACF if you already have computed it.')
    general.add_argument('-s','--save',type=str,default=None,
                         help='if -s option is used it will save a numpy file of the ACF, default not to save')
    general.add_argument('--oned', type=str, default=None,
                         help='give --oned f or --oned t if you want to compute the one dimensional ACF in frequency or time, respectively. For frequency it will sum between -t and -T to determine the spectrum. For time it will sum across the entire frequency band. Default is to do a 2D ACF.')
    general.add_argument('--fscrunch',type=int, default=1,
                         help='Give fscrunch factor to downsample the archive')
    general.add_argument('--tscrunch',type=int, default=1,
                         help='Give tscrunch factor to downsample the archive')
    general.add_argument('-D', '--distance', type=float, default=None,
                         help='Distance to the source in Mpc')
    general.add_argument('-S', '--sefd', type=float, default=20/1.54,
                         help='system equivalent flux density [default: 20/1.54 for Eff]')
    parser.add_argument('--db', type=str, default=None,
                        help='Path to pickle file that contains the pandas data frame.')
    parser.add_argument('--src', type=str, default='any',
                        help='Used when writing to a pandas data frame. In case parameters ' +
                        'for the same burst are determined in different ways (e.g. different files)' +
                        'this parameter can be used to distinguish them.')
    parser.add_argument('--compare_to', type=str, default=None,
                        help='Used when writing to a pandas data frame and things are to be' +
                        'compared to a different --src entry. If that other entry exists, all' +
                        'fitting ranges and similar from that one will be used.')
    parser.add_argument('--newfit', action='store_true',
                        help='If set will prompt the user to identify the' +
                        'regions in time for the ACF. Only used if --df is used.')
    parser.add_argument('--show_acf_plots', action='store_true',
                        help='If set will show diagnotic ACF plots.')
    return parser.parse_args()


if __name__ == "__main__":
    args = options()

    # load in archive
    archive, extent, tsamp, T0 = load_archive(args.arfile, dm=args.dm,
                                              pscrunch=True, extent=True,
                                              fscrunch=args.fscrunch,
                                              tscrunch=args.tscrunch)
    print(archive.shape)
    t_min, t_max, f_min, f_max = extent
    if f_min > f_max:
        f = f_min
        f_min = f_max
        f_max = f
    fsamp = (f_max - f_min) / archive.shape[0]  # in MHz
    BW = f_max - f_min
    # tsamp comes in seconds

    # if there is a database with central pixels we use the initial values from there
    if args.db is not None:
        src = args.src
        compare_to = args.compare_to
        if not os.path.exists(args.db):
            data, ncolumns = init_pandas_df()
        else:
            data = pd.read_pickle(args.db)
            ncolumns = data.columns.shape[0]
        try:
            exp, dish, scan = (os.path.basename(args.arfile)).split('_')[0:3]
        except:
            raise ValueError('File name does not agree with assumed convention.')
        df = data[(data.experiment == exp) &
                  (data.dish == dish) &
                  (data.scan == scan) &
                  (data.src == src)]
        nexisting = 0
        if not df.empty:
            print(f'Found entry for mode {src} in {exp}_{dish}_{scan}.')
            central_bins = np.array(df.t0_ms.values / (tsamp * 1e3))
            nexisting = len(central_bins)
            central_freqs = np.array(df.f0_MHz.values)
            #print(central_bins)
            ranges = np.array([[b, e] for b, e in zip(df.begin_acf_range.values,
                                                      df.end_acf_range.values)]) / (tsamp * 1e3)
            #print(ranges)
            # in case this a new entry we certainly want find new ranges.
            for Range in ranges:
                b, e = Range
                if b == e == 0:
                    args.newfit = True
            off_range_t = np.array([df.off_range_t0.values[0],
                                    df.off_range_t1.values[0]]) // args.tscrunch

        else:
            print(f'No entry for mode {src} in {exp}_{dish}_{scan} '+
                  f'in {args.db}. Generating a new one.')
            args.newfit = True
        if compare_to is not None:
            df_comp = data[(data.experiment == exp) &
                           (data.dish == dish) &
                           (data.scan == scan) &
                           (data.src == compare_to)]
            if not df_comp.empty:
                args.newfit = False
                print(f'Found entry for mode {compare_to} in {exp}_{dish}_{scan}.')
                central_bins = np.array(df_comp.t0_ms.values / (tsamp * 1e3))
                central_freqs = np.array(df_comp.f0_MHz.values)
                nexisting = len(central_bins)
                #print(central_bins)
                ranges = np.array([[b, e] for b, e in zip(df_comp.begin_acf_range.values,
                                                          df_comp.end_acf_range.values)]) / (tsamp * 1e3)
                off_range_t = np.array([df_comp.off_range_t0.values[0],
                                        df_comp.off_range_t1.values[0]]) // args.tscrunch
                for component in range(len(central_bins)):
                    # if entry exists update it, if it doesn't exist add a new one.
                    if df[(df.component == component)].empty:
                        df_list = [0] * ncolumns
                        df_list[:5] = exp, dish, scan, component, src
                        newrow, _ = init_pandas_df(vals=[df_list])
                        data = pd.concat([data, newrow])
                    b, e = ranges[component]
                    peak = central_bins[component]
                    while (b > peak) or (e < peak):
                        subbursts = identify_bursts(np.mean(archive, axis=0),
                                                    f"The peak (={peak}) is outside the chosen range for component {component}. Pick the peak via holding down i and clicking left on the peak")

                        central_bins[component] = subbursts.peak_times[0]
                    #print(ranges[component])
                    data.loc[(data.experiment == exp) &
                             (data.dish == dish) &
                             (data.scan == scan) &
                             (data.src == src) &
                             (data.component == component),
                             't0_ms'] = central_bins[component] * (tsamp * 1e3)
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             'begin_acf_range'] = ranges[component][0] * (tsamp * 1e3)
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             'end_acf_range'] = ranges[component][1] * (tsamp * 1e3)
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             'off_range_t0'] = off_range_t[0] * args.tscrunch
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             'off_range_t1'] = off_range_t[1] * args.tscrunch
            else:
                print(f'Found no entry for {compare_to} in {exp}_{dish}_{scan}.')
                args.newfit = True
        if args.newfit:
            overview = identify_bursts_2d(archive, tsamp, 'Overview',
                                          block=False)
            print(f'Please identify the edges and peaks of subcomponents for the ACF by '+
                  'left clicking in the plot. Peaks via i + left click.')
            subbursts = identify_bursts(np.mean(archive, axis=0),
                                        "Identify edges and peaks of subcomponents for the ACF.")
            ranges = np.array(subbursts.ranges).reshape(-1, 2)
            central_bins = np.array(subbursts.peak_times)
            subbursts = identify_bursts(np.mean(archive, axis=0),
                                        "Please choose the off-pulse region (in time).")
            off_range_t = np.array(subbursts.ranges, dtype=np.int)

            if (len(central_bins) == 0) & (len(ranges) == 0):
                raise ValueError(f'Seems you did not pick anything?')
            if not len(central_bins) == len(ranges):
                raise ValueError(f'The selected number of ranges does not agree with the number of components.')

            # in case we have entries from an earlier attempt where we used more
            # components than this time, we'll delete the ones that are 'too many'
            if nexisting > len(central_bins):
                data_exp = data.loc[(data.experiment == exp) &
                                    (data.dish == dish) &
                                    (data.src == src) &
                                    (data.scan == scan)]
                data = data.loc[(data.experiment != exp) |
                                (data.dish != dish) |
                                (data.src != src) |
                                (data.scan != scan)]
                for component in range(len(central_bins), nexisting):
                    data_exp = data_exp.loc[(data_exp.component != component)]
                data = pd.concat([data, data_exp])

            for component in range(len(central_bins)):
                # if entry exists update it, if it doesn't exist add a new one.
                if df[(df.component == component)].empty:
                    df_list = [0] * ncolumns
                    df_list[:5] = exp, dish, scan, component, src
                    newrow, _ = init_pandas_df(vals=[df_list])
                    data = pd.concat([data, newrow])
                b, e = ranges[component]
                peak = central_bins[component]
                while (b > peak) or (e < peak):
                    subbursts = identify_bursts(np.mean(archive, axis=0),
                                                f"The peak (={peak}) is outside the chosen range for component {component}. Pick the peak via holding down i and clicking left on the peak")

                    central_bins[component] = subbursts.peak_times[0]
                #print(ranges[component])
                data.loc[(data.experiment == exp) &
                         (data.dish == dish) &
                         (data.scan == scan) &
                         (data.src == src) &
                         (data.component == component),
                         't0_ms'] = central_bins[component] * (tsamp * 1e3)
                data.loc[(data.experiment == exp) & (data.dish == dish) &
                         (data.src == src) &
                         (data.scan == scan) & (data.component == component),
                         'begin_acf_range'] = ranges[component][0] * (tsamp * 1e3)
                data.loc[(data.experiment == exp) & (data.dish == dish) &
                         (data.src == src) &
                         (data.scan == scan) & (data.component == component),
                         'end_acf_range'] = ranges[component][1] * (tsamp * 1e3)
                data.loc[(data.experiment == exp) & (data.dish == dish) &
                         (data.src == src) &
                         (data.scan == scan) & (data.component == component),
                         'off_range_t0'] = off_range_t[0] * args.tscrunch
                data.loc[(data.experiment == exp) & (data.dish == dish) &
                         (data.src == src) &
                         (data.scan == scan) & (data.component == component),
                         'off_range_t1'] = off_range_t[1] * args.tscrunch
    else:
        begintime = args.begintime
        if args.endtime != 0:
            if args.endtime <= args.begintime:
                raise ValueError('Please give an end time after the begin time, i.e. -T > -t')
            else:
                endtime = args.endtime
        else:
            subbursts = identify_bursts(np.mean(archive,axis=0),
                                        "Identify edges of subcomponents for the ACF.")
            ranges = np.array(subbursts.ranges).reshape(-1, 2)

    archive_orig = archive.copy()
    for component, Range in enumerate(ranges):
        archive = archive_orig.copy()
        # at this point ranges are floats
        begintime, endtime = Range
        begintime = round(begintime)
        endtime = round(endtime)
        archive = archive[:, begintime:endtime]


        # compute the ACF
        if args.loadnpy is None and args.oned is None:
            ACF=autocorr_2D(archive)
        elif args.loadnpy is not None:
            ACF = np.load(args.loadnpy)
        elif args.oned == 't' and args.loadnpy is None:
            profile = np.mean(archive, axis=0)
            ACF = autocorr_fft(profile)  # , len(profile), v=None)
            ACF = np.concatenate((ACF[::-1], ACF[1:]))
        elif args.oned == 'f' and args.loadnpy is None:
            spec = np.mean(archive,axis=1)
            mask = np.ones_like(spec)
            mask[np.where(spec == 0)[0]] = 0
            ACF = autocorr_fft(spec)  # ,len(spec),v=mask)
            ACF = np.concatenate((ACF[::-1],ACF[1:]))

        ACF/=np.max(ACF)
        #mask the zero lag spike
        ACFmasked = np.ma.masked_where(ACF==np.max(ACF),ACF)

        fit_mask = ACFmasked.mask
        fit_mask = fit_mask.astype(np.float)

        ACFoff = autocorr_fft(np.mean(archive_orig[:,0:endtime-begintime],axis=0))   #8:10],axis=1))
        ACFoff = np.concatenate((ACFoff[::-1],ACFoff[1:]))
        ACFoff/=np.max(ACFoff)
        ACFoffmasked = np.ma.masked_where(ACFoff==np.max(ACFoff),ACFoff)

        print(len(ACFoff))
        print(len(ACFmasked))

        if args.save is not None:
            np.save(args.save,ACF)

        if args.oned is not None:
            #define axes
            time_archive = np.arange(1,archive.shape[1],1)*tsamp*1000 #ms
            times = np.concatenate((-time_archive[::-1],np.concatenate(([0],time_archive))))
            freq_res = (extent[3]-extent[2])/archive.shape[0] #MHz
            freq_archive = np.arange(0,archive.shape[0]-1,1)*freq_res
            freqs = np.concatenate((-freq_archive[::-1],np.concatenate(([0],freq_archive))))

            gmodel = Model(lorentz)
            if args.oned  == 'f':
                ACF_for_fit = ACFmasked[int(len(ACF)/2.-(10/freq_res)):int(len(ACF)/2.+(10/freq_res))]
                x_for_fit = freqs[int(len(ACF)/2.-(10/freq_res)):int(len(ACF)/2.+(10/freq_res))]
                x=freqs
                print("*** Lorentzian fit to freq ACF ***")
                result = gmodel.fit(ACF_for_fit, x=x_for_fit, gamma=0.00005, y0=0.01, c=0)
            if args.oned == 't':
                ACF_for_fit = ACFmasked[int(len(ACF)/2.-(0.05/(tsamp*1000))):int(len(ACF)/2.+(0.05/(tsamp*1000)))]
                x_for_fit = times[int(len(ACF)/2.-(0.05/(tsamp*1000))):int(len(ACF)/2.+(0.05/(tsamp*1000)))]
                x=times
                print("*** Lorentzian fit to time ACF ***")
                result = gmodel.fit(ACF_for_fit, x=x_for_fit, gamma=0.005, y0=1, c=0)

            print(result.fit_report())
            plt.plot(x, ACFmasked,'k',alpha=0.8)
            plt.plot(x,lorentz(x,result.params['gamma'],result.params['y0'],result.params['c']),color='purple')
            plt.show()

            plt.plot(x,ACFmasked-lorentz(x,result.params['gamma'],result.params['y0'],result.params['c']),'k',alpha=0.8)
            plt.show()


        if args.oned is None:
            ACFtime = np.sum(ACF,axis=0)
            ACFfreq = np.sum(ACF,axis=1)
            ACFtime = np.ma.masked_where(ACFtime==np.max(ACFtime),ACFtime)
            ACFfreq = np.ma.masked_where(ACFfreq==np.max(ACFfreq),ACFfreq)

            time_archive = np.arange(1,archive.shape[1],1)*tsamp*1000 #ms
            times = np.concatenate((-time_archive[::-1],np.concatenate(([0],time_archive))))
            freq_res = (extent[3]-extent[2])/archive.shape[0] #MHz
            freq_archive = np.arange(0,archive.shape[0]-1,1)*freq_res
            freqs = np.concatenate((-freq_archive[::-1],np.concatenate(([0],freq_archive))))

            #1D Gaussian fitting to ACFtime and ACF freq
            time_ACF = freq_ACF = True
            P0t = [1,0,np.max(times)]
            P0f = [1,0,np.max(freqs)]
            try:
                poptt, pcovt = curve_fit(gaus, times, ACFtime, p0=P0t)
            except:
                print('Gaussian fit to time ACF failed.')
                time_ACF = False
                poptt = P0t
            try:
                poptf, pcovf = curve_fit(gaus, freqs, ACFfreq, p0=P0f)
            except:
                print('Gaussian fit to freq ACF failed.')
                freq_ACF = False
                poptf = P0f
            if (not time_ACF) or (not freq_ACF):
                args.show_acf_plots = True

            if args.show_acf_plots:
                fit=False
                while fit==False:
                    fig,ax=plt.subplots(2)
                    ax[0].plot(times, ACFtime)
                    ax[0].plot(times,gaus(times,*poptt))
                    ax[0].set_ylabel('ACF time')
                    ax[1].plot(freqs, ACFfreq)
                    ax[1].plot(freqs,gaus(freqs,*poptf))
                    ax[1].set_ylabel('ACF freq')
                    plt.show()

                    print(poptt,pcovt)
                    answer = input("Are you happy with this fit? (y/n): ")
                    if answer == 'n':
                        guess = input("Give an initial guess (time_amp,time_mean,time_sigma,freq_amp,freq_mean,freq_sigma) ")
                        poptt, pcovt = curve_fit(gaus, times, ACFtime, p0=[guess[0],guess[1],guess[2]])
                        poptf, pcovf = curve_fit(gaus, freqs, ACFfreq, p0=[guess[3],guess[4],guess[5]])
                    if answer == 'y':
                        fit=True

            #2D Gaussian fitting
            timesh, freqs_m = np.meshgrid(times, freqs)
            timesh = timesh.astype('float64')
            freqs_m = freqs_m.astype('float64')

            """
            data = twoD_Gaussian((timesh,freqs_m), 3, 0, 0, poptt[2],poptf[2],0)
            data=data.reshape(len(freqs),len(times))
            plt.imshow(data,aspect='auto',origin='lower')
            plt.show()
            exit()
            """

            params = Parameters()
            params.add('amplitude', value=1)
            params.add('xo',value=0,vary=False)
            params.add('yo',value=0,vary=False)
            params.add('sigma_x',value=poptt[2],min=poptt[2]-0.2*poptt[2], max=poptt[2]+0.2*poptt[2])
            params.add('sigma_y',value=poptf[2],min=poptf[2]-0.2*poptf[2], max=poptf[2]+0.2*poptf[2])
            params.add('theta',value=0)



            #twoD_Gaussian_fit(params,x_data_tuple,data)
            out = minimize(twoD_Gaussian_fit, params, kws={"x_data_tuple": (timesh,freqs_m), "data": ACFmasked})
            print("*** Gaussian fit to 2D ACF ***")
            print(fit_report(out))
            if args.db is not None:
                sigma_x = out.params['sigma_x'].value  # time width in ms
                sigma_y = out.params['sigma_y'].value  # freq width in MHz
                sigma_x_err = out.params['sigma_x'].stderr
                sigma_y_err = out.params['sigma_y'].stderr
                hpbw_f = 2.355 * np.abs(sigma_y)
                if (compare_to is None) or args.newfit:
                    if hpbw_f*2 < BW*0.75:
                        # chose central freq, else use middle
                        subbursts = identify_bursts_2d(archive_orig, tsamp,
                                                       f"Pick the central frequency bin for component {component}.")
                        central_freq = subbursts.peak_freqs[0] * fsamp + f_min
                        use_full_BW = False
                    else:
                        central_freq = f_min + BW / 2.
                        use_full_BW = True
                else:
                    central_freq = central_freqs[component]
                    use_full_BW = False if (hpbw_f*2 < BW*0.75) else True

                print(f'Using central frequency of {central_freq} MHz')
                for col, val in zip(['twidth_sigma_ms', 'fwidth_sigma_MHz',
                                     'twidtherr_ms', 'fwidtherr_MHz', 'f0_MHz'],
                                    [sigma_x, sigma_y, sigma_x_err,
                                     sigma_y_err, central_freq]):
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             col] = val

            #popt, pcov = curve_fit(twoD_Gaussian, (timesh, freqs_m), ACFmasked.flatten().astype('float64') , p0=guess, sigma=std_array)
            #data_fitted = twoD_Gaussian((timesh, freqs_m), *popt)
            data_fitted = twoD_Gaussian((timesh, freqs_m), out.params['amplitude'],out.params['xo'],out.params['yo'],out.params['sigma_x'],out.params['sigma_y'],out.params['theta'])
            data_fitted = data_fitted.reshape(len(freqs),len(times))

            #residuals


            ACFtimeresid = ACFtime-np.sum(data_fitted,axis=0)
            ACFfreqresid = ACFfreq-np.sum(data_fitted,axis=1)
            ACFt_for_fit = ACFtimeresid[int(len(ACFtimeresid)/2.-(0.01/tsamp)):int(len(ACFtimeresid)/2.+(0.01/tsamp))]
            ACFf_for_fit = ACFfreqresid[int(len(ACFfreqresid)/2.-(20/freq_res)):int(len(ACFfreqresid)/2.+(20/freq_res))]
            time_for_fit = times[int(len(ACFtimeresid)/2.-(0.01/tsamp)):int(len(ACFtimeresid)/2.+(0.01/tsamp))]
            freq_for_fit = freqs[int(len(ACFfreqresid)/2.-(20/freq_res)):int(len(ACFfreqresid)/2.+(20/freq_res))]

            #popt_lorentz_time, pcov_lorentz_time = curve_fit(lorentz, time_for_fit, ACFt_for_fit,p0=[0.0005,1,0])
            #popt_lorentz_freq, pcov_lorentz_freq = curve_fit(lorentz, freq_for_fit, ACFf_for_fit,p0=[1,1,0])

            gmodel = Model(lorentz)
            result_time = gmodel.fit(ACFt_for_fit, x=time_for_fit, gamma=0.001, y0=1, c=0)
            result_freq = gmodel.fit(ACFf_for_fit, x=freq_for_fit, gamma=1, y0=1, c=0)
            print("*** Lorentzian fit to time ACF ***")
            print(result_time.fit_report())
            print("*** Lorentzian fit to freq ACF ***")
            print(result_freq.fit_report())

            if args.db is not None:
                scint_bw = result_freq.params['gamma'].value  # time width in ms
                scint_bw_err = result_freq.params['gamma'].stderr  # time width in ms
                for col, val in zip(['scint_bw_MHz', 'scint_bw_err_MHz'],
                                    [scint_bw, scint_bw_err]):
                    data.loc[(data.experiment == exp) & (data.dish == dish) &
                             (data.src == src) &
                             (data.scan == scan) & (data.component == component),
                             col] = val

            if args.save!=None:
                np.save("fit"+str(args.save),np.array([times,freqs,data_fitted,lorentz(freqs,result_freq.params['gamma'],result_freq.params['y0'],result_freq.params['c'])]))
            #plot
            if args.show_acf_plots:
                fig = plt.figure(figsize=(8, 8))
                rows=3
                cols=3
                widths = [3, 1,1]
                heights = [1,1,3]
                gs = gridspec.GridSpec(ncols=cols, nrows=rows,width_ratios=widths, height_ratios=heights, wspace=0.0, hspace=0.0)

                cmap = plt.cm.gist_yarg

                ax1 = fig.add_subplot(gs[0,0]) # Time ACF
                ax1.plot(times,ACFtime,color='k')
                plt.setp(ax1.get_xticklabels(), visible=False)
                plt.setp(ax1.get_yticklabels(), visible=False)
                ax1.set_xlim(times[0],times[-1])
                ax1.plot(times,np.sum(data_fitted,axis=0),color='purple')

                ax2 = fig.add_subplot(gs[1,0],sharex=ax1) # Time ACF residuals
                ax2.plot(times,ACFtimeresid,color='k')
                plt.setp(ax2.get_yticklabels(), visible=False)
                plt.setp(ax2.get_xticklabels(), visible=False)
                ax2.set_xlim(times[0],times[-1])
                ax2.plot(times, lorentz(times,result_time.params['gamma'],result_time.params['y0'],result_time.params['c']),color='orange')

                ax3 = fig.add_subplot(gs[2,0],sharex=ax2) # 2D ACF
                T,F=np.meshgrid(times, freqs)
                ax3.imshow(ACFmasked,aspect='auto',interpolation='nearest',origin='lower',cmap=cmap,extent=(times[0],times[-1],freqs[0],freqs[-1]))
                ax3.contour(T,F,data_fitted,4, colors='r', linewidths=.5)
                ax3.set_ylabel('Freq lag [MHz]')
                ax3.set_xlabel('Time lag [ms]')

                ax4 = fig.add_subplot(gs[2,1],sharey=ax3) #Freq ACF residuals
                ax4.plot(ACFfreqresid,freqs,color='k')
                plt.setp(ax4.get_yticklabels(), visible=False)
                plt.setp(ax4.get_xticklabels(), visible=False)
                ax4.set_ylim(freqs[0],freqs[-1])
                ax4.plot(lorentz(freqs,result_freq.params['gamma'],result_freq.params['y0'],result_freq.params['c']),freqs,color='orange')

                ax5 = fig.add_subplot(gs[2,2],sharey=ax4) #Freq ACF
                ax5.plot(ACFfreq,freqs,color='k')
                plt.setp(ax5.get_yticklabels(), visible=False)
                plt.setp(ax5.get_xticklabels(), visible=False)
                ax5.set_ylim(freqs[0],freqs[-1])
                ax5.plot(np.sum(data_fitted,axis=1),freqs,color='purple')

                plt.show()

        if use_full_BW:
            beginchan, endchan = 0, archive_orig.shape[0]-1
        else:
            beginchan = round((central_freq - f_min - hpbw_f) / fsamp)
            endchan = round((central_freq - f_min + hpbw_f) / fsamp)
            beginchan = 0 if beginchan < 0 else beginchan
            endchan = archive_orig.shape[0]-1 if endchan >= archive_orig.shape[0] else endchan
        if off_range_t[0] == off_range_t[1] == 0:
            print(f"Please choose the off-pulse region (in time).")
            subbursts = identify_bursts(np.mean(archive_orig, axis=0),
                                        "Please choose the off-pulse region (in time).")
            off_range_t = np.array(subbursts.ranges, dtype=np.int)
        print(beginchan, endchan, begintime, endtime)
        archive_chans_full_t = archive_orig[beginchan:endchan, :].copy()
        archive = archive[beginchan:endchan, :]

        profile = np.mean(archive, axis=0)
        profile_full = np.mean(archive_chans_full_t, axis=0)

        profile_off = np.mean(archive_chans_full_t[:, off_range_t[0]:off_range_t[1]],
                              axis=0)

        #convert to S/N
        profile -= np.mean(profile_off)
        profile_full -= np.mean(profile_off)
        profile_off -= np.mean(profile_off)
        profile /= np.std(profile_off)
        profile_full /= np.std(profile_off)

        rows = 2
        cols = 1
        gs = gridspec.GridSpec(rows, cols, wspace=0., hspace=0.,
                               height_ratios=[0.5,]*(rows-1) + [2,],
                               width_ratios=[5,] + [1,]*(cols-1))

        pltrange = [int(begintime-30//(tsamp*1e3)), int(endtime+50//(tsamp*1e3))]
        prof = plt.subplot(gs[0])
        prof.plot(profile_full)
        prof.set_ylabel('S/N')
        #plt.axvline(off_range_t[0], color='red', label='off-region')
        #plt.axvline(off_range_t[1], color='red')
        prof.axvline(begintime, color='blue', label='on-region')
        prof.axvline(endtime, color='blue')
        prof.legend()
        prof.set_xlim(pltrange)
        #plt.show()
        #plt.clf()
        vmin = archive_orig.std() * -1
        vmax = archive_orig.std() * 6
        dyn = plt.subplot(gs[1], sharex=prof)
        dyn.imshow(archive_orig, aspect='auto', interpolation='Gaussian',
                   origin='lower', vmin=vmin, vmax=vmax)
        dyn.axhline(beginchan, color='green')
        dyn.axhline(endchan, color='green')
        #plt.axvline((t_min + off_range_t[0]*tsamp) * 1000, color='red', label='off-region')
        #plt.axvline((t_min + off_range_t[1]*tsamp) * 1000, color='red')
        dyn.axvline(begintime, color='red', label='on-region')
        dyn.axvline(endtime, color='red')
        dyn.set_xlim(prof.get_xlim())
        plt.show()
        #plt.clf()
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
