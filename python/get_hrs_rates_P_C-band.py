import pandas as pd
import numpy as np


# function which can take overlap into account between scans
# taken from Mark Snelders
def merge(times):
    """ Function which reduces the overlap between a list of tuples
    Example [(1, 3), (2, 4), (7, 8)] ---> [(1, 4), (7, 8)]
    Use: unique_times = list(merge(times)) """
    saved = list(times[0])
    for st, en in sorted([sorted(t) for t in times]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


def count_bursts(ids, threshhold, spc, sfxc):
    '''
    ids: list of strings with burst IDs
    threshold: threshold above which to return bursts (in Jyms)
    spc: pandas data frame that has spc-corrected fluences
    spxc: pandas data frame that has sfxc-measured fluences (used only for stockert bursts)
    '''
    ids.sort()
    bursts = {}
    counter = 0
    id_nrs = [' '] # initialise with one element
    for ID in ids:
        id_nr = ID[:3]
        # The spc ones
        tot_fluence = np.sum(spc[(spc.id == ID)].fluence_jyms.values)
        if 'st' in ID:
            tot_fluence = np.sum(sfxc_sub[(sfxc_sub.id == ID)].fluence_jyms.values)
        if tot_fluence >= threshhold:
            bursts[ID] = round(tot_fluence, 2)
            counter += 1
            if id_nr == id_nrs[-1]:
                counter -= 1
        id_nrs.append(id_nr)
    return counter, bursts


def expected_counts(hrs_with_detections, counts, hrs_without_detections):
    return hrs_without_detections / hrs_with_detections * counts


df = pd.read_pickle("../dbs/burst_info.pickle")
obs_df = pd.read_pickle("../dbs/df_r67_reduced.pkl")

sfxc = df[(df.src == "sfxc")] # contains the ToAs
spc = df[(df.src == "spc")]   # for our 'best' fluences

# FAST window
mjdstart = 59305
mjdstop = 59363

# 2nd activity window
mjdstart = 59602
mjdstop = 59642

# to get the observing hours
obs_df_sub = obs_df[(obs_df.mjd_start >= mjdstart) &
                    (obs_df.mjd_end < mjdstop)]

# and the associated bursts
sfxc_sub = sfxc[(sfxc.toa_bary_tdb_inf_freq >= mjdstart) &
                (sfxc.toa_bary_tdb_inf_freq < mjdstop)]
ids = sfxc_sub.id.unique()

# P-band hrs
freq_max = 400  #in MHz
p_hrs = np.sum(obs_df_sub[(obs_df_sub.cent_f < freq_max)].rec_time_h.values)
print(f'P-band hours: {p_hrs:.2f} hrs')

# C-band hrs
freq_min = 2300  #in MHz
c_hrs = np.sum(obs_df_sub[(obs_df_sub.cent_f > freq_min)].rec_time_h.values)
print(f'C-band hours: {c_hrs:.2f} hrs')

# L-band hrs, here we risk overlap with another dish
freq_min = 1000  #in MHz
freq_max = 2000  #in MHz
l_obs = obs_df_sub[(obs_df_sub.cent_f > freq_min) &
                   (obs_df_sub.cent_f < freq_max)]
l_obs = l_obs.sort_values(by=['mjd_start'])
start_mjds = l_obs.mjd_start
end_mjds = l_obs.mjd_end
l_hrs = 0
times = [(x, y) for x, y in zip(l_obs.mjd_start.values, l_obs.mjd_end.values)]
times_red = list(merge(times))
for i in range(len(times_red)):
    s, e = times_red[i][0], times_red[i][1]
    l_hrs += (e-s)
l_hrs *= 24
print(f'L-band hours: {l_hrs:.2f} hrs (non-overlapping)')


# 15-sigma detection threshold, in Jyms
thresh_P = 91
thresh_C = 5

# asuming a flat spectrum
counts, pbursts = count_bursts(ids, thresh_P, spc, sfxc)
expected = expected_counts(l_hrs, counts, p_hrs)
#print(pbursts)
print(f'\nFound a total of {counts} bursts above threshold of {thresh_P} Jyms, assuming a flat spectrum!')
print(f'Thus we would have expected {expected:.1f} bursts at P-band.')

counts, cbursts = count_bursts(ids, thresh_C, spc, sfxc)
expected = expected_counts(l_hrs, counts, c_hrs)

#print(cbursts)
print(f'\nFound a total of {counts} bursts above threshold of {thresh_C} Jyms, assuming a flat spectrum!')
print(f'Thus we would have expected {expected:.1f} bursts at C-band.')

# if spectral index is different from 0, say -1.5, then the threshold of
# detections at L-band is lower/higher. We assume the same width.
# spectral index of emission
alpha = -1.5
L = 1424 # in MHz, typical central freq at L-band (O8)
P = 332  # in MHz, typical central freq at P-band (Wb)
C = 4678 # in MHz, typical central freq at C-band (Tr)

thresh_P = thresh_P * (L/P)**(alpha)
thresh_C = thresh_C * (L/C)**(alpha)

counts, pbursts = count_bursts(ids, thresh_P, spc, sfxc)
expected = expected_counts(l_hrs, counts, p_hrs)
#print(pbursts)
print(f'\nFound a total of {counts} bursts above threshold of {thresh_P:.1f} Jyms, assuming alpha = {alpha}!')
print(f'Thus we would have expected {expected:.1f} bursts at P-band.')

counts, cbursts = count_bursts(ids, thresh_C, spc, sfxc)
expected = expected_counts(l_hrs, counts, c_hrs)

#print(cbursts)
print(f'\nFound a total of {counts} bursts above threshold of {thresh_C:.1f} Jyms, assuming alpha = {alpha}!')
print(f'Thus we would have expected {expected:.1f} bursts at C-band.')
