import argparse
import numpy as np
import psrchive
import pandas as pd
from load_file import load_filterbank
from presto import filterbank
from interactive_sub_burst_identifier_1d import identify_bursts
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy import coordinates as coord
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time import Time


def init_src_station(ra, dec, location='coe', xyz=None):
    stations = {'coe': [0., 0., 0.],
                'o8': [3370965.8787, 711466.1978, 5349664.2006],
                'tr': [3638558.5100, 1221969.7200, 5077036.7600],
                'wb': [3828750.6969, 442589.2176, 5064921.5700],
                'st': [4031510.647, 475159.114, 4903597.840]
                }
    if xyz is not None:
        assert len(xyz) == 3
        x, y, z = xyz
    else:
        if location not in stations.keys():
            raise ValueError(f'Location {location} not implemented.' +
                             f'Choices are {stations.keys()}')
        else:
            x, y, z = stations[location]
    source = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    station = EarthLocation.from_geocentric(x=x * u.m,
                                            y=y * u.m,
                                            z=z * u.m)
    return source, station


def multi_gauss(params):
    amp, x0, sigma = params[0]
    model = models.Gaussian1D(amp, x0, sigma)
    if len(params) > 1:
        for param in params[1:]:
            amp, x0, sigma = param
            model += models.Gaussian1D(amp, x0, sigma)
    return model


def fit_gauss(model, data, times):
    fit_LM = fitting.LevMarLSQFitter()
    fit = fit_LM(model, times, data)
    return fit, fit_LM


def float2str(filcoord):
    neg = ''
    if filcoord < 0:
        filcoord *= -1
        neg = '-'
    ci = int(filcoord)
    cf = str(filcoord - ci).split('.')[1]
    ci = f'{ci:06d}'
    h, m, s = [f'{ci[i:i+2]}' for i in range(0, len(ci), 2)]
    return f'{neg}{h}:{m}:{s}.{cf}'


def get_fil_coords(filfile):
    ra = filfile.header['src_raj']
    dec = filfile.header['src_dej']
    ra = float2str(ra)
    dec = float2str(dec)
    return ra, dec


def calc_time_delay(freq, dm, dmconst):
    """
    Input:
        freq (MHz), float
        dm (pc/cc), float
        dmconst (GHz^2 cm^3 pc^-1 ms), float
    Returns:
        The time difference between infinite frequency and freq
        in SECONDS"""
    return dmconst * 10**6. * dm * freq**-2. / 1000.


def loc2bary(mjd, dm, ref_freq, location, source, dm_const=4.14880568679703):
    localTOA = Time(mjd, format='mjd', scale='utc', location=location)
    delay_to_inf = calc_time_delay(ref_freq, dm, dm_const) * u.s  # in seconds
    delay_to_bary = localTOA.light_travel_time(source, 'barycentric')  # in TDB seconds, to be added
    baryTOA = localTOA - delay_to_inf + delay_to_bary
    return baryTOA.tdb.value


def get_flags_T0_from_archive(archive_name):
    print("Loading in archive file {}".format(archive_name))
    archive = psrchive.Archive_load(archive_name)
    w = archive.get_weights().squeeze()
    w[w>0] = 1
    # the channel numbering in archives is flipped compared to filterbanks
    w = np.flip(w)
    tsamp = archive.get_first_Integration().get_duration()/archive.get_nbin()
    T0 = archive.get_first_Integration().get_start_time().in_days()  # comes as MJD
    return w.astype(np.bool), T0, tsamp


def window(ds,window,tsamp):
    """
    chop out a window around the burst (peak of the profile)
    """
    profile = np.mean(ds,axis=0)
    begin=np.argmax(profile)-int(window/(1000*tsamp))
    end=np.argmax(profile)+int(window/(1000*tsamp))
    burstds=ds[:,begin:end]
    return burstds,begin


def downsamp(ds, tdown=1, fdown=1):
    if fdown != 1:
        ds = ds.reshape(ds.shape[0]//fdown, fdown, ds.shape[-1]).sum(axis=1)
    if tdown != 1:
        ds = ds.reshape(ds.shape[0], ds.shape[-1]/tdown, tdown).sum(axis=2)
    return ds


def convert_SN(burst_prof, off_prof):
    burst_prof -= np.mean(off_prof)
    off_prof -= np.mean(off_prof)
    burst_prof /= np.std(off_prof)
    return burst_prof


def get_clean_dyn(dynspec, dynspec_off, flag_mask=None):
    # remember, in the mask True indicates invalid, i.e. a masked entry
    dyn = np.zeros_like(dynspec)
    #removing bandpass
    for fr in range(dynspec.shape[0]):
        dyn[fr, :] = convert_SN(dynspec[fr, :], dynspec_off[fr, :])
    if flag_mask is not None:
        dyn = np.ma.masked_array(dyn, mask=~flag_mask)
    return dyn


def options():
    parser = argparse.ArgumentParser()
    general = parser.add_argument_group('General info about the data.')
    general.add_argument('filterbank', type=str, default=None,
                         help='Filterbank file to perform the Gaussian fits on.')
    general.add_argument('-i', '--ID', type=str, required=True,
                         help='Unique burst identifier in database.')
    general.add_argument('-r', '--half_range', type=float, default=0.5,
                         help='half range to extract around burst. In units of seconds.')
    general.add_argument('--arfile', type=str, default=None,
                         help='Archive file that contains the same burst and that was used to compute the fluences.')
    general.add_argument('-d','--dm',type=float, default=0,
                         help='DM to dedisperse the data to.')
    general.add_argument('-t', '--burst_time', type=float, default=None,
                         help='Time of burst into the filterbank in seconds.')
    general.add_argument('--fscrunch',type=int, default=1,
                         help='Give fscrunch factor to downsample the archive')
    general.add_argument('--tscrunch',type=int, default=1,
                         help='Give tscrunch factor to downsample the archive')
    parser.add_argument('--db', type=str, default=None, required=True,
                        help='Path to pickle file that contains the pandas data frame.')
    parser.add_argument('--src', type=str, default='sfxc',
                        help='Used when writing to a pandas data frame. In case parameters ' +
                        'for the same burst are determined in different ways (e.g. different files)' +
                        'this parameter can be used to distinguish them.')
    parser.add_argument('--location', type=str, default='coe',
                        help='Location of TOA, i.e. location for which times in the filterbank' +
                        'are true, e.g. Tr, O8, St. Default is COE (i.e. geocentric)')
    parser.add_argument('--ref_freq', type=float, default=None,
                        help='Reference frequency in MHz to which the data are/were dedispersed. ' +
                        'Only used to compute the barycentric TOA at infinite frequency.' +
                        'Useful in case the data are already (coherently) dedispersed.' +
                        'If unused will take the top of the band as per the header of the filterbank.')
    parser.add_argument('--show_dynspec', action='store_true',
                        help='If set will show a diagnostic plot of the on-pulse dynamic spectrum and profile.')
    return parser.parse_args()


if __name__ == "__main__":
    args = options()
    db = pd.read_pickle(args.db)
    ID = args.ID
    filname = args.filterbank
    src = args.src
    archive = args.arfile
    half_range = args.half_range
    dm = args.dm
    burst_time = args.burst_time
    tscrunch = args.tscrunch
    location = args.location
    ref_freq_MHz = args.ref_freq
    show_dynspec = args.show_dynspec
    df = db[(db.id == ID) & (db.src == src)]
    if df.empty:
        raise ValueError(f'No entry for {ID} and src {src}.')
    peaks = np.array(df.t0_ms.values) / 1e3  # seconds
    twidths = np.array(df.gauss2d_twidth_ms.values) / 1e3  # seconds
    dm_4_bary = dm
    if dm == 0.0:
        dm_4_bary = float(df.dm.values[0])

    if archive is not None:
        flags, T0_ar, tsamp_ar = get_flags_T0_from_archive(archive)
        if burst_time is None:
            burst_time = T0_ar + peaks.mean()/86400.  # MJD
    if burst_time is None:
        # TODO: this is dodgy..
        raise ValueError(f'burst_time is required if no archive file is supplied.')

    # need to estimate how far into the filterbank our burst lives
    fil = filterbank.FilterbankFile(filname)
    T0_fil = fil.header['tstart']
    dur = fil.nspec * fil.header['tsamp']
    ra = float2str(fil.header['src_raj'])
    dec = float2str(fil.header['src_dej'])
    fil.close()
    source, station = init_src_station(ra, dec, location)
    secs_in = (burst_time - T0_fil) * 86400  # seconds
    if args.burst_time is not None:
        secs_in = args.burst_time
    try:
        assert secs_in > 0
    except:
        print(f'Time into filterbank is negative: secs_in = {secs_in}')
        print(f'Time from archive: {burst_time}, '+
              f'Filterbank start: {T0_fil}.')
        print(f'Will load the whole thing.')
        half_range = dur / 2
        secs_in = half_range

    print(f'Seeking {secs_in} into the file, extracting +/- {half_range}s around that time.')
    arr, extent, tsamp_fil, fsamp_fil, begbin = load_filterbank(filname, dm=dm,
                                                                burst_time=secs_in,
                                                                half_range=half_range,
                                                                tscrunch=tscrunch)
    print(f'First bin is {begbin} which corresponds to {tsamp_fil*begbin} seconds into the file.')
    f0 = extent[2]
    if ref_freq_MHz is None:
        ref_freq_MHz = extent[3]

    if flags is None:
        flags = np.ones_like(arr, dtype=np.bool)
    print(flags.dtype, flags.sum(), flags.shape)
    arr = np.ma.masked_array(arr)
    arr[flags == 0] = np.ma.masked
    subbursts = identify_bursts(arr.mean(axis=0),
                                "Please choose the off-pulse region (in time).")
    off_range = np.array(subbursts.ranges, dtype = np.int)
    arr = get_clean_dyn(arr,
                        arr[:, off_range[0]:off_range[1]])
    #print(arr.shape)
    profile_full = arr.mean(axis=0)
    if show_dynspec:
        #print(profile_full.shape)
        rows = 2
        cols = 1
        gs = gridspec.GridSpec(rows, cols, wspace=0., hspace=0.,
                               height_ratios=[0.5,]*(rows-1) + [2,],
                               width_ratios=[5,] + [1,]*(cols-1))

        prof = plt.subplot(gs[0])
        prof.plot(profile_full)
        prof.set_ylabel('S/N')
        prof.legend()
        vmin = arr.std() * -1
        vmax = arr.std() * 6
        dyn = plt.subplot(gs[1], sharex=prof)
        dyn.imshow(arr, aspect='auto', interpolation='Gaussian',
                   origin='lower', vmin=vmin, vmax=vmax)
        plt.show()
        plt.clf()

    subbursts = identify_bursts(profile_full,
                                f"Identify the fit-range, i.e. the on-time.")
    on_range = np.array(subbursts.ranges, dtype=np.int)
    answer = 'n'
    profile_org = profile_full.copy()
    tscrunch_tot = 1
    while answer == 'n':
        moveon = False
        while not moveon:
            subbursts = identify_bursts(profile_full, pltrange=on_range,
                                        title=f"Identify the peaks of components ({len(peaks)}) and the widths.")
            widths = np.array(subbursts.ranges).reshape(-1, 2)
            widths = np.array([(comp[1]-comp[0])/4 for comp in widths])
            central_bins = np.array(subbursts.peak_times)
            central_amps = np.array(subbursts.peak_amps)
            if len(central_bins) == len(peaks):
                moveon = True
            else:
                print(f'Selected number of peaks does not agree with what is in the data base. Expected {len(peaks)} peaks.')
                moveon = False


        # twidths come in seconds per component -- use as guess for width
        # central_bins are from filterbank in units of bins
        # on_range also comes in units of bins
        # and everthing is relative to begbin from the chopping
        # since most things are in bin-space, fit in that, then convert to times.
        amps = central_amps
        begintime = on_range[0] #* tsamp_fil
        endtime = on_range[1] #* tsamp_fil
        central_bins -= begintime
        profile_fit = profile_full[begintime:endtime]
        guess = np.array([[amp, x0, sigma] for amp, x0, sigma in zip(amps, central_bins, widths)])
        model = multi_gauss(guess)
        bins = np.arange(len(profile_fit))
        fit, fit_LM = fit_gauss(model, profile_fit, bins)
        plt.plot(profile_fit)
        plt.plot(fit(bins))
        #plt.xlim(pltrange)
        plt.show()
        # ask if you're happy with the fit
        # if not shall we either downsample
        # or shall we pick new peaks
        # return the fit, the downsampled profile
        answer = input("Are you happy with this fit? (y/n): ")
        if answer == 'n':
            tscrunch = input("Shall we downsample? By which factor?(power of 2!) Or shall we pick new peaks? (p) Or abort (a)?")
            if tscrunch == 'a':
                quit(1)
            elif tscrunch == 'p':
                tscrunch = 1
            elif int(tscrunch) % 2 == 0:
                tscrunch = int(tscrunch)
            else:
                raise ValueError(f'Input invalid. Can be a, p or power of 2.')
            if tscrunch > 1:
                rest = len(profile_full) % tscrunch
                profile_full = profile_full[0:len(profile_full)-rest]
                profile_full = profile_full.reshape(-1, tscrunch).mean(axis=1)
                on_range //= tscrunch
                tscrunch_tot *= tscrunch

    central_bins = (fit.parameters.reshape(-1, 3)[:, 1] + begintime) * tscrunch_tot + begbin
    sigmas = (fit.parameters.reshape(-1, 3)[:, 2]) * tscrunch_tot * tsamp_fil * 1e3
    TOAs = central_bins * tsamp_fil / 86400. + T0_fil
    for i, TOA in enumerate(TOAs):
        print(f'TOA for component {i} = {TOA}')
        db.loc[(db.id == ID) &
               (db.src == src) &
               (db.component == i),
               'toa_mjd_utc'] = TOA
        db.loc[(db.id == ID) &
               (db.src == src) &
               (db.component == i),
               'toa_loc'] = location
        db.loc[(db.id == ID) &
               (db.src == src) &
               (db.component == i),
               'ref_freq_MHz'] = ref_freq_MHz
        db.loc[(db.id == ID) &
               (db.src == src) &
               (db.component == i),
               'toa_bary_tdb_inf_freq'] = loc2bary(TOA, dm_4_bary, ref_freq_MHz, station, source)
        db.loc[(db.id == ID) &
               (db.src == src) &
               (db.component == i),
               'gauss2d_twidth_ms'] = sigmas[i]
    db.to_pickle(args.db)
