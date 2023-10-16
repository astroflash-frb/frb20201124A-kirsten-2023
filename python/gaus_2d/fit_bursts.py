"""
Kenzie Nimmo 2020
Version to take archive file and perform 2D fit

"""
import sys
import numpy as np
import os
import optparse
import pandas as pd
import matplotlib.pyplot as plt
from load_file import load_archive, load_filterbank
import burst_2d_Gaussian as fitter

from scipy.stats import chisquare
from interactive_sub_burst_identifier_2d import identify_bursts as identify_bursts_2d
from interactive_sub_burst_identifier_2d import ds
from interactive_sub_burst_identifier_1d import identify_bursts as identify_bursts_1d


def save_fit_params(pandasframe, exp, dish, scan, params, errors,
                    start_bin, tsamp, T0):
    '''
    Meant to either append fit params to a common data base. If that pandas data frame
    does not yet exist we'll creat it.
    '''
    columns = ['experiment', 'dish',
               'scan', 'component',
               't0_mjd', 'f0_MHz',
               'terr_ms', 'ferr_MHz',
               'gauss2d_twidth_ms', 'gauss2d_fwidth_MHz',
               'begin_acf_range', 'end_acf_range',
               'twidth_sigma_ms', 'fwidth_sigma_MHz',
               'twidtherr_ms', 'fwidtherr_MHz',
               'off_range_t0', 'off_range_t1',
               'fluence_jyms', 'peak_snr', 'peak_flux_jy',
               'scint_bw_MHz', 'scint_bw_err_MHz',
               'energy_iso_erg_per_hz', 'spectral_lum_erg_per_s_per_hz']
    if not os.path.exists(pandasframe):
        print(f'Initiating new pandas pickle file at {pandasframe}')
        data = pd.DataFrame(columns=columns)
    else:
        data = pd.read_pickle(pandasframe)
        print('loaded')
    # check if the entry already exists, if so keep what's in it but replace
    # with new values
    for component, param in enumerate(params):
        # t0 from the fit comes in units of millisecond from the start of the
        # fitting window, tsamp comes in units of seconds.
        t0_mjd = param[1] + start_bin * (tsamp * 1e3)
        f0_mhz = param[2]
        t0err = errors[component, 1]
        f0err = errors[component, 2]
        t_width = param[3]
        f_width = param[4]

        if not data[(data.experiment == exp) & (data.dish == dish) &
                    (data.scan == scan) & (data.component == component)].empty:

            print(f'Updating entry for: {exp}_{dish}_{scan} component {component}.')
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     't0_mjd'] = t0_mjd
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     'f0_MHz'] = f0_mhz
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     'terr_ms'] = t0err
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     'ferr_MHz'] = f0err
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     'gauss2d_twidth_ms'] = t_width
            data.loc[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan) & (data.component == component),
                     'gauss2d_fwidth_MHz'] = f_width
        else:
            print(f'Adding new entry for: {exp}_{dish}_{scan} component {component}.')
            df_list = [[exp, dish, scan, component, t0_mjd, f0_mhz, t0err, f0err,
                        t_width, f_width, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
            df = pd.DataFrame(df_list, columns=columns)
            data = pd.concat([data, df])
    # in case we started anew and now have fewer components than in an
    # earlier run, we delete those other components
    ncomponents = len(params)
    nexisting = data[(data.experiment == exp) & (data.dish == dish) &
                     (data.scan == scan)].component.count()
    if nexisting > ncomponents:
        data_exp = data.loc[(data.experiment == exp) & (data.dish == dish) &
                            (data.scan == scan)]
        data = data.loc[(data.experiment != exp) & (data.dish != dish) &
                        (data.scan != scan)]
        for component in range(ncomponents, nexisting):
            data_exp = data_exp.loc[(data_exp.component != component)]
        data = pd.concat([data, data_exp])
    data.to_pickle(pandasframe)
    return

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [options] infile', description="2D Gaussian fit "
                                   "to FRB data. Input an archive file.")
    parser.add_option('-t', '--tavg', dest='tavg', type='int', default=1,
                      help="If -t option is used, time averaging is applied using the factor "
                           "given after -t.")
    parser.add_option('-s', '--subb', type='int', default=1,
                      help="If -s option is used, subbanding is applied using the factor "
                           "given after -s.")
    parser.add_option('-f', '--tfit', dest='tfit', type='float',
                      help="Give the time before and after the burst in ms that should be used "
                      "for fitting.")
    parser.add_option('--ptavg', type='int', default=1,
                      help="As -t but only used for plotting.")
    parser.add_option('--psubb', type='int', default=1,
                      help="As -s but only used for plotting.")
    parser.add_option('-d','--dm', type='float', default=None,
                      help="Give the DM to correct the archive to.")
    parser.add_option('-g', '--pguess', action='store_true',
                      help="Plot the guessed functions")
    parser.add_option('-y', '--yes', action='store_true',
                      help="Automatically give y for all questions.")
    parser.add_option('--skip', action='store_true',
                      help="Skip second round of fitting with drifting Gaussian.")
    parser.add_option('--filterbank', action='store_true',
                      help="If you want to read in a filterbank instead of an archive file use this option.")
    parser.add_option('--fullpol', action='store_true',
                      help="If the filterbank file contains full polarisation give this option in addition to --filterbank")
    parser.add_option('--db', type=str, default=None,
                      help='Path to pickle file that contains the pandas data frame. If the file ' +
                      'does not exist it will be created.')
    options, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)
    elif len(args) != 1:
        sys.stderr.write("Only one input file must be provided!\n")
    else:
        options.infile = args[-1]

    tavg = options.tavg
    subb = options.subb
    ptavg = options.ptavg
    psubb = options.psubb

    #read in the data
    basename = options.infile

    if options.filterbank is None:
        print("2D Gaussian fit of archive %s" % (basename))
    else:
        print("2D Gaussian fit of filterbank %s" % (basename))

    if options.dm!=None:
        dm = options.dm
    else:
        dm = 0

    if options.filterbank is None:
        waterfall, extent, tsamp, T0 = load_archive(basename, dm=dm, tscrunch=tavg,
                                                    fscrunch=subb, extent=True,
                                                    pscrunch=True,
                                                    remove_baseline=True)
        print(T0, extent, tsamp)
    else:
        if options.fullpol is not None:
            print("Not dedispersing the filterbank file since you specified full pol")
            waterfall, extent, tsamp = load_filterbank(basename, dm=None,
                                                     fullpol=True)
        else:
            waterfall, extent, tsamp = load_filterbank(basename, dm=dm,
                                                     fullpol=False)
    #waterfall[:, 1330:1500] = 0
    profile = np.mean(waterfall, axis=0)

    print("Pick the range to be shown in the 2D-plot for the burst picking.")
    subbursts = identify_bursts_1d(profile, "Pick the range to be shown in the 2D-plot for the burst picking.")
    begin_samp, end_samp = np.array(subbursts.ranges, dtype=np.int)

    # remove_baseline=True in load_archive has calibrated the bandpass
    waterfall = waterfall[:, begin_samp:end_samp]
    subbursts = identify_bursts_2d(waterfall, tsamp)
    time_guesses = np.array(subbursts.peak_times)
    freq_guesses = np.array(subbursts.peak_freqs, dtype=np.float32)
    amp_guesses = np.array(subbursts.peak_amps)

    # Get time before and after the burst to be used for fitting.
    if options.tfit:
        tfit = options.tfit
    else:
        tfit = 5

    fit_mask = np.zeros(waterfall.shape, dtype=np.bool)
    fit_start = int(np.min(subbursts.peak_times) - tfit*1e-3/tsamp)
    if fit_start < 0:
        fit_start = 0

    fit_end = int(np.max(subbursts.peak_times) + tfit*1e-3/tsamp)
    if fit_end > waterfall.shape[1]:
        fit_end = waterfall.shape[1]-1
    fit_mask[:, fit_start:fit_end] = True
    fit_mask[np.where(waterfall == 0)] = False
    fit_mask = fit_mask.astype(np.float)

    freqs = np.linspace(extent[2], extent[3], waterfall.shape[0])
    fsamp_orig = freqs[1]-freqs[0]
    freqs_orig = freqs.copy()

    if subb != 1:
        waterfall = ds(waterfall, factor=subb, axis=0)
        fit_mask = ds(fit_mask, factor=subb, axis=0)
        freqs = ds(freqs, factor=subb, axis=0)
        fsamp=fsamp_orig*subb

    # Get the times at the pixel centers in ms.
    times = (np.arange(waterfall.shape[1]) * tsamp + tsamp/2) * 1e3
    #print(times.shape, fit_start, fit_end)
    start_stop = times[[fit_start, fit_end-1]]

    time_guesses *= tsamp*1e3
    #freq_guesses= np.array(map(float, freq_guesses))
    freq_guesses *= fsamp_orig
    freq_guesses += freqs_orig[0]

    n_sbs = len(time_guesses)#waterfall.shape[0]
    freq_std_guess = [50.] * n_sbs  # comes in MHz. n_sbs * [int(512 / subb / 4.)]
    t_std_guess = [1.0] * n_sbs  # comes in ms
    #amp_guesses = np.sqrt(tavg*subb) * amp_guesses #/ np.sqrt((~waterfall.mask[:,0]).sum())

    use_standard_2D_gaussian = True
    fix_angle = False
    while True:
        if use_standard_2D_gaussian:
            print("Using standard Gaussian model")
            model = fitter.gen_Gauss2D_model(time_guesses, amp_guesses,
                                             f0=freq_guesses,
                                             bw=freq_std_guess, dt=t_std_guess, verbose=True)
            if fix_angle:
                angles = model.param_names[5::6]
                for a in angles:
                    model.fixed[a] = True
        else:
            # Use the custom drifting Gaussian model.
            print("Using drifting Gaussian model")
            model = fitter.drifting_2DGaussian(
                amplitude=amp_guesses[0], t_mean=time_guesses[0], f_mean=freq_guesses[0],
                t_stddev=t_std_guess[0], f_stddev=freq_std_guess[0], drift=0.)
            for a, tg, fg, tsg, fsg in zip(amp_guesses[1:], time_guesses[1:],
                                           freq_guesses[1:], t_std_guess[1:],
                                           freq_std_guess[1:]):
                model += fitter.drifting_2DGaussian(amplitude=a, t_mean=tg, f_mean=fg,
                                                    t_stddev=tsg, f_stddev=fsg, drift=0.)
            if fix_angle:
                angles = model.param_names[5::6]
                for a in angles:
                    model.fixed[a] = True

        # Plot the guess (only once)
        if options.pguess and use_standard_2D_gaussian:
            low_res_waterfaller = ds(ds(waterfall, psubb), ptavg, axis=1)
            fitter.plot_burst_windows(ds(times, ptavg), ds(freqs, psubb), low_res_waterfaller,
                                      model, ncontour=8, res_plot=True, vlines=start_stop)

        bestfit, fitLM = fitter.fit_Gauss2D_model(waterfall, times, freqs, model,
                                                  weights=fit_mask)
        bestfit_params, bestfit_errors, corr_fig = fitter.report_Gauss_parameters(bestfit,
                                                                                  fitLM,
                                                                                  verbose=True)

        timesh, freqs_m = np.meshgrid(times, freqs)
        fit_mask = fit_mask.astype(np.bool)
        timesh, freqs_m = timesh[fit_mask], freqs_m[fit_mask]
        chisq, pvalue = chisquare(waterfall[fit_mask], f_exp=bestfit(timesh, freqs_m),
                                  ddof=6*n_sbs, axis=None)
        print("Chi^2 and p-value:", chisq, pvalue)
        # Plot downsampled and subbanded as given in options
        low_res_waterfaller = ds(ds(waterfall, psubb), ptavg, axis=1)
        fig, res_fig = fitter.plot_burst_windows(ds(times, ptavg), ds(freqs, psubb),
                                                 low_res_waterfaller, bestfit, ncontour=8,
                                                 res_plot=True, vlines=start_stop)  # diagnostic plots
        if use_standard_2D_gaussian:
            corr_fig.savefig('%s_correlation.pdf'%basename)
            fig.savefig('%s_fit.pdf'%basename)
            res_fig.savefig('%s_fit_residuals.pdf'%basename)

        if options.yes:
            answer = 'y'
        else:
            plt.show()
            print("Note the begin bin of the window is "+str(begin_samp))

            answer = input("Are you happy with the fit y/n/skip/fix/fb?")

        if (answer == 'y') and use_standard_2D_gaussian:
            use_standard_2D_gaussian = False
            if options.skip:
                if options.db is not None:
                    try:
                        exp, dish, scan = (os.path.basename(basename)).split('_')[0:3]
                    except:
                        raise ValueError('File name does not agree with assumed convention.')
                    save_fit_params(options.db, exp, dish, scan,
                                    bestfit_params, bestfit_errors, begin_samp,
                                    tsamp, T0)
                break
        elif answer == 'n':
            print("Try and use different guesses this time.")
            use_standard_2D_gaussian = True
        elif (answer == 'skip'):
            if options.db is not None:
                try:
                    exp, dish, scan = (os.path.basename(basename)).split('_')[0:3]
                except:
                    raise ValueError('File name does not agree with assumed convention.')
                save_fit_params(options.db, exp, dish, scan,
                                bestfit_params, bestfit_errors, begin_samp,
                                tsamp, T0)
                break
        elif answer == 'fix':
            fix_angle = not fix_angle
            if fix_angle:
                print("The angle is fixed now.")
            else:
                print("The angle is a fittable parameter again.")
