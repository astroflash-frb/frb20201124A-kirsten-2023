import pandas as pd
import numpy as np


df = pd.read_pickle("../dbs/burst_info.pickle")
obs_df = pd.read_pickle("../dbs/df_r67_reduced.pkl")

sfxc = df[(df.src == "sfxc")] # contains the ToAs

mjdstarts = [59305, 59602]
mjdstops = [59363, 59642]
freq_min = 1000  #in MHz
freq_max = 2000  #in MHz
stations = ['o8', 'wb', 'tr', 'stk', 'st']

for mjdstart, mjdstop in zip(mjdstarts, mjdstops):
    print(f'Epoch 1: {mjdstart}--{mjdstop}')
    # to get the observing hours at L-band
    obs_df_sub = obs_df[(obs_df.mjd_start >= mjdstart) &
                        (obs_df.mjd_end < mjdstop) &
                        (obs_df.cent_f > freq_min) &
                        (obs_df.cent_f < freq_max)]

    # and the associated bursts
    sfxc_sub = sfxc[(sfxc.toa_bary_tdb_inf_freq >= mjdstart) &
                    (sfxc.toa_bary_tdb_inf_freq < mjdstop)]
    for station in stations:
        hrs = round(obs_df_sub[(obs_df_sub.Telescope) == station].rec_time_h.values.sum(), 0)
        nbursts = len(sfxc_sub[(sfxc_sub.dish) == station].id.unique())

        print(f'For {station} we have {nbursts} burst in {hrs} hours.')
