import pandas as pd


def init_pandas_df(vals=None):
    columns = ['experiment',
               'dish',
               'scan',
               'component',
               'src',
               't0_ms',
               'f0_MHz',
               'toa_loc',
               'terr_ms',
               'ferr_MHz',
               'gauss2d_twidth_ms',
               'gauss2d_fwidth_MHz',
               'begin_acf_range',
               'end_acf_range',
               'twidth_sigma_ms',
               'fwidth_sigma_MHz',
               'twidtherr_ms',
               'fwidtherr_MHz',
               'off_range_t0',
               'off_range_t1',
               'fluence_jyms',
               'peak_snr',
               'peak_flux_jy',
               'fluence_tot_jyms',
               'scint_bw_MHz',
               'scint_bw_err_MHz',
               'energy_iso_erg_per_hz',
               'spectral_lum_erg_per_s_per_hz',
               'id',
               'sefd_jy',
               'distance_Mpc',
               'dm',
               'toa_mjd_utc',
               'ref_freq_MHz',
               'toa_bary_tdb_inf_freq']

    if vals is None:
        return pd.DataFrame(columns=columns), len(columns)
    else:
        assert len(vals) == 1
        assert len(vals[0]) == len(columns)
        return pd.DataFrame(vals, columns=columns), len(columns)
