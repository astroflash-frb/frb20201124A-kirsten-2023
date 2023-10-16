import psrchive
import numpy as np
from presto import filterbank


def load_archive(archive_name, rm=0, dm=0, tscrunch=None, fscrunch=None,
                 remove_baseline=True, extent=False, model=False,
                 pscrunch=False, cent=False):
    print("Loading in archive file {}".format(archive_name))
    archive = psrchive.Archive_load(archive_name)
    if cent is True:
        archive.centre()

    archive.tscrunch()
    if pscrunch is True:
        archive.pscrunch()

    ardm = archive.get_dispersion_measure()
    ardd = archive.get_dedispersed()
    arfc = archive.get_faraday_corrected()
    arrm = archive.get_rotation_measure()

    if ardd is True:
        print("Archive file is already dedispersed to a DM of {} pc/cc".format(ardm))
        if dm != 0:
            dm = dm - ardm

    if remove_baseline is True:
        #remove the unpolarised background -- note this step can cause Stokes I < Linear
        archive.remove_baseline()

    if (tscrunch is not None) and (tscrunch != 1):
        archive.bscrunch(tscrunch)

    if (fscrunch is not None) and (fscrunch != 1):
        archive.fscrunch(fscrunch)

    if dm != 0:
        if ardd is True and dm == 0:
            pass
        if ardd is True and ardm != dm:
            print("Dedispersing the data to a DM of {} pc/cc".format(dm+ardm))
        else:
            print("Dedispersing the data to a DM of {} pc/cc".format(dm))
        archive.set_dispersion_measure(dm)
        archive.set_dedispersed(False)
        archive.dedisperse()

    if rm != 0:
        if arfc is True:
            print("Faraday derotating the data to a RM of {} rad/m^2".format(rm+arrm))
        else:
            print("Faraday derotating the data to a RM of {} rad/m^2".format(rm))
        archive.set_rotation_measure(rm)
        archive.set_faraday_corrected(False)
        archive.defaraday()

    ds = archive.get_data().squeeze()
    w = archive.get_weights().squeeze()
    print(ds.shape, w.shape)
    if model is False:
        if len(ds.shape) == 3:
            weights = np.ones_like(ds)
            weights[w == 0, :, :] = 0
            ds = np.ma.masked_where(weights == 0, ds)
            ds = np.flip(ds, axis=1)
        else:
            weights = np.ones_like(ds)
            weights[w == 0, :] = 0
            ds = np.ma.masked_where(weights == 0, ds)
            ds = np.flip(ds, axis=0)

    if model is True:
        ds = ds

    tsamp = archive.get_first_Integration().get_duration()/archive.get_nbin()

    if extent is True:
        extent = [0, archive.get_first_Integration().get_duration()*1000,
                  archive.get_centre_frequency()+archive.get_bandwidth()/2.,
                  archive.get_centre_frequency()-archive.get_bandwidth()/2.]
        T0 = archive.get_first_Integration().get_start_time().in_days()
        return ds, extent, tsamp, T0
    else:
        return ds


def loaddata(filename, t_burst, DM=0, maskfile=None, fullpol=False, window=10):
    """
    Reads in the filterbank file into a numpy array, applies the mask
    and outputs the burst dynamic spectrum as a numpy array.
    Inputs:
    - filename: Name of filterbank filename with full path
    - t_burst: time of the burst in seconds into the filterbank file
    - maskfile: text file containing a list of frequency channels to zap (get it from pazi command)
    - fullpol: True/False. full polarisation filterbanks are not currently supported in this analysis. Stokes I only.
    - window: window in [ms] around burst to extract for analysis (default +-10ms)
    Outputs:
    - Stokes I dynamic spectrum
    - Off burst dynamic spectrum
    - time resolution in seconds
    - frequency resolution in MHz
    """
    if not fullpol:
        ds, dsoff, extent, tsamp, begbin = load_filterbank(filename,
                                                           dm=DM,
                                                           fullpol=False,
                                                           burst_time=t_burst)
        StokesI_ds = np.zeros_like(ds)
        StokesI_off = np.zeros_like(dsoff)
        #removing bandpass
        for fr in range(ds.shape[0]):
            StokesI_ds[fr, :] = convert_SN(ds[fr, :], dsoff[fr, :])
            StokesI_off[fr, :] = convert_SN(dsoff[fr, :], dsoff[fr, :])

        # frequency resolution
        freqres = (extent[3]-extent[2])/ds.shape[0]
        # frequency array
        frequencies = np.linspace(extent[2], extent[3], ds.shape[0])

        if maskfile!=None:
            maskchans=np.loadtxt(maskfile, dtype='int')
            maskchans = [StokesI_ds.shape[0]-1-x for x in maskchans]
            StokesI_ds[maskchans, :] = 0
            StokesI_off[maskchans, :] = 0

    else:
        raise ValueError('Full pol filterbank data is not currently supported.')

    return StokesI_ds, StokesI_off, tsamp, freqres, begbin, frequencies


def convert_SN(burst_prof, off_prof):
    burst_prof -= np.mean(off_prof)
    off_prof -= np.mean(off_prof)
    burst_prof /= np.std(off_prof)
    return burst_prof


def load_filterbank(filterbank_name, dm=0.0, fullpol=False, burst_time=None,
                    half_range=None, tscrunch=1):
    '''
    Shamelessly stolen and modified from Kenzie's bbpipe.py
    '''
    fil = filterbank.FilterbankFile(filterbank_name)
    tsamp = fil.header['tsamp']
    # consider how much delay the DM would cause
    tdel = np.abs(4.149377593360996e-3*dm*((1./np.min(fil.frequencies/1000.)**2)
                                           - (1./np.max(fil.frequencies/1000.)**2)))  # seconds

    if burst_time is not None:
        # if you know where the burst is in seconds within the filterbank file,
        # chop out only necessary chunk of data
        burst_bin = int(burst_time/tsamp)
        if half_range is not None:
            nspec = int(half_range / tsamp)
            begbin = burst_bin - nspec
            endbin = burst_bin + nspec
            begbin = 0 if begbin < 0 else begbin
            endbin = fil.nspec-10 if endbin >= fil.nspec else endbin
            spec = fil.get_spectra(begbin, endbin)
        else:
            if burst_bin < fil.nspec//2:
                off = fil.get_spectra(fil.nspec//2 + 100, fil.nspec-10)
            else:
                off = fil.get_spectra(0, (fil.nspec//2)-100)

            # burst is at 50ms unless burst is within 50ms of the start of the file
            if (burst_bin-int(50e-3/tsamp)+int(2*tdel/tsamp) < fil.nspec) & \
               (burst_bin-int(50e-3/tsamp) >= 0):
                spec = fil.get_spectra(burst_bin-int(50e-3/tsamp),
                                       int(2*tdel/tsamp))
                begbin = burst_bin-int(50e-3/tsamp)
            elif burst_bin-int(50e-3/tsamp) < 0:
                spec = fil.get_spectra(0, int(2*tdel/tsamp))
                begbin = 0
            else:
                dur = fil.nspec - (burst_bin-int(50e-3/tsamp))
                spec = fil.get_spectra(burst_bin-int(50e-3/tsamp), dur)
                begbin = burst_bin-int(50e-3/tsamp)

    else:
        spec = fil.get_spectra(0, fil.nspec)
        begbin = 0

    if burst_time is not None and fullpol:
        raise ValueError("Definining the burst time in ful pol data is not currently supported")
    if dm is not None and fullpol:
        raise ValueError("If filterbank contains full polarisation information, dedispersion won't work properly")
    if dm is not None:
        spec.dedisperse(dm)
        if (burst_time is not None) and (half_range is None):
            off.dedisperse(dm)
            offarr = off.data

    arr = spec.data
    if (burst_time is not None) and (dm > 0.0):
        #chop off the end where the DM delay creates an angled edge
        arr = arr[:, :-int(tdel/tsamp)]
    if fullpol:
        arr = arr.reshape(fil.header['nchans'], -1, 4)

    if fil.header['foff'] < 0:
        # this means the band is flipped
        arr = np.flip(arr, axis=0)
        foff = fil.header['foff']*-1
        if (burst_time is not None) and (half_range is None):
            offarr = np.flip(offarr, axis=0)
    else:
        foff = fil.header['foff']

    #header information
    tsamp = fil.header['tsamp']
    begintime = 0
    endtime = arr.shape[1]*tsamp
    fch_top = fil.header['fch1']
    nchans = fil.header['nchans']
    fch_bottom = fch_top+foff-(nchans*foff)

    if tscrunch > 1:
        rest = arr.shape[1] % tscrunch
        arr = arr[:, 0:arr.shape[1]-rest].reshape(arr.shape[0], -1, tscrunch).mean(axis=-1)
        endtime -= rest * tsamp
        tsamp *= tscrunch
        begbin //= tscrunch
    extent = (begintime, endtime, fch_bottom, fch_top)

    if burst_time is None:
        return arr, extent, tsamp
    if (burst_time is not None) and (half_range is None):
        return arr, offarr, extent, tsamp, begbin
    else:
        return arr, extent, tsamp, foff, begbin
