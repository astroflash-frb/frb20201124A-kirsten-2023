import pandas as pd
import numpy as np
from collections import OrderedDict as od


# actual ontime or total width? Width gives a lover limit on the luminosity.
def speclum(fluence, ontime, distance=453, z=0.098):
    '''
    fluence assumed in Jy~ms
    ontime assumed in ms
    distance assumed in Mpc
    '''
    #convert Jy ms to J s; and ms to s
    fluence_jys = fluence*1e-3
    ontime /= 1e3
    #convert Mpc to cm
    distance_lum_cm = 3.086e24*distance
    energy_iso= fluence_jys*4*np.pi*(distance_lum_cm**2)*1e-23 / (1+z)**(2)
    lum_spec = energy_iso/ontime
    return lum_spec

## Fluence as sum of components
## t-Width as total width, i.e. start(component_0) - end(component_n)
## spectral_lum computed from sum(fluence), done on the fly.
##

# Columns
#           ID  Station TOA       peak_SNR  Fluence             nComponents  t-Width Spectral_lum                 BW   scint_BW
template0 = '{0} & {1} & {2:.9f} & {3:.1f}   & ${4:.1f}\pm{5:.1f}$ & {6}   &    {7:.2f}   & ${8:.2f}\pm{9:.2f}$  & {10} & {11} \\\ \n'
# for bursts without baseband data we add a * to the MJD as a different DM constant was used
template1 = '{0} & {1} & {2:.9f}* & {3:.1f}   & ${4:.1f}\pm{5:.1f}$ & {6}   &    {7:.2f}   & ${8:.2f}\pm{9:.2f}$ & {10} & {11} \\\ \n'

header = '\\begin{table*}\n' +\
         '\\caption{\\label{tab:burst_properties}Burst properties.}\n' +\
         '\\footnotesize \n' +\
         '\\begin{tabular}{lccccccccc}' +\
         '\\hline\n' +\
         '\\hline\n'

columns = od()
columns[1] = ['Burst ID & ',
              ' & ']
columns[2] = ['Station & ',
              ' & ']
columns[3] = [r'{TOA$\mathrm{^{a}}$} & ',
              r'{[MJD]} & ']
columns[4] = [r'{Peak S/N$\mathrm{^{b}}$} & ',
              ' & ']
columns[5] = [r'{Fluence$\mathrm{^{c}}$} & ',
              r'{[Jy~ms]} & ']
columns[6] = ['Number of & ',
              'components & ']
columns[7] = [r'{Width$\mathrm{^{d}}$} & ',
              r'{[ms]} & ']
columns[8] = [r'{Spectral Luminosity$\mathrm{^{e}}$} & ',
              r'{[$\mathrm{10^{32}\,erg\,s^{-1}\,Hz^{-1}}$]} & ']
columns[9] = [r'BW$\mathrm{^{f}}$ & ',
             r'{[MHz]} & ']
columns[10] = [r'{$\nu_{s}\mathrm{^{g}}$}' + '\\\\ \n ',
              r'{[MHz]}' + '\\\\ \n']

heading = ''
# first row of heading
for key, value in columns.items():
    heading += value[0]
# second row of heading
for key, value in columns.items():
    heading += value[1]

footer = '\\hline\n' +\
    '\multicolumn{9}{l}{$\mathrm{^{a}}$ Time of arrival at the solar system barycenter at infinite frequency in TDB (using a DM of $410.8$~\dmunit, a dispersion measure} \\\ \n' +\
    '\multicolumn{9}{l}{\hspace{0.25cm}constant of $1/0.000241$ $\mathrm{GHz^2~cm^3~pc^{-1}~\mu s}$ and as (J2000) position RA = $05\!:\!08\!:\!03.5$, Dec = $+26\!:\!03\!:\!37.8$; for times} \\\ \n' +\
    '\multicolumn{9}{l}{\hspace{0.25cm}marked with a *, a DM constant of $4.14880568679703$ $\mathrm{GHz^2\,cm^3\,pc^{-1}\,ms}$ was used).} \\\ \n' +\
    '\multicolumn{9}{l}{\hspace{0.25cm}For multi-component busts, the TOA is defined as the middle between the peak of the first and the last component.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{b}}$ The peak S/N of the brightest component.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{c}}$ Computed as the sum over the measured fluence of each component. We assume a conservative error of $20$\% for all bursts, dominated by the} \\\ \n' +\
    '\multicolumn{9}{l}{\hspace{0.25cm}uncertainty of the SEFD.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{d}}$ Manually determined time span between start of first and end of last component.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{e}}$ Computed using $D_L=453~\mathrm{Mpc}$, $z=0.098$ and the listed width.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{f}}$ Bandwidth used for computing the fluence. This is often the full available observing bandwidth.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{g}}$ Weighted average over the measured scintillation bandwidth per component.} \\\ \n' +\
    '\multicolumn{9}{l}{$\mathrm{^{h}}$ No measurement was possible due to low S/N.} \\\ \n' +\
    '\\hline\n' +\
    '\\end{tabular}\n' +\
    '\\normalsize\n' +\
    '\\end{table*}\n\n\n'

df = pd.read_pickle('../dbs/burst_info.pickle')
# We want to list only the SPC'ed fluences. That's fine for most bursts but not for the Stockert ones
# and for B04-o8. Stockert has only src=sfxc entries and B04-o8 has both src=scale and src=sfxc, where the
# latter is a copy of the former for plotting reasons (we plot the SFXC generated versions of all bursts.
# To have them all we concatenate the SPC-frames for most with the SFXC-frames for Stockert and the SCALE-frame
# for B04-o8, then we
# - order them in time
# - compute the middle between first and last component for the overall TOA
# - sum up the fluences
# - compute the spectral luminosity from there
df_mjd = df[(df.src == 'sfxc')] # this contains the TOAs and is used as standalone df

# these three will be concatenated to get all burst fluences
df_spc = df[(df.src == 'spc')]  # this is for the 'correct' fluences of o8, wb, tr
df_st = df[(df.dish == 'st')]   # stockert has no spc-entry
df_B04 = df[(df.id == 'B04-o8') & (df.src == 'scale')] # this burst has no 'spc' entry because there are no baseband data
df_B31 = df[(df.id == 'B31-o8') & (df.src == 'scale')] # this burst has no 'spc' entry because scaling did not work

df = pd.concat([df_spc, df_st, df_B04, df_B31])
df_mjd.sort_values(by='toa_bary_tdb_inf_freq', ignore_index=True, inplace=True)
ids = df_mjd.id.unique()

bw98 = ['B01-st', 'B06-st', 'B10-st', 'B31-st', 'B32-st', 'B33-st', 'B34-st', 'B35-st',
        'B37-st', 'B38-st', 'B39-st', 'B40-st', 'B43-st', 'B44-st', 'B45-st', 'B46-st']
bw256 = ['B08-tr', 'B22-o8']
bw512 = ['B25-o8', 'B26-o8', 'B27-o8', 'B31-o8', 'B32-o8']
scintbws = []
fluence_100 = 0
fluence_500 = 0
table = ''
for ID in ids:
    bw = 128
    if ID in bw98:
        bw = 98
    if ID in bw256:
        bw = 256
    if ID in bw512:
        bw = 512
    burst = df[(df.id == ID)]
    mjd_frame = df_mjd[(df_mjd.id == ID)]
    toas = mjd_frame.toa_bary_tdb_inf_freq.values
    toa = (toas.min() + toas.max()) / 2.
    ref_freq = mjd_frame.ref_freq_MHz.values[0]
    cent_f = ref_freq - bw/2.
    dish = (burst.dish.unique()[0]).capitalize()
    ncomp = len(burst.component.values)
    peakSNR = burst.peak_snr.values.max()
    fluence = burst.fluence_jyms.values.sum()
    signal_bw = int(abs(burst.fwidth_sigma_MHz.values.max()) * 2.355 * 2)
    if signal_bw > bw * 0.75:
        signal_bw = bw
    if fluence >= 100.0:
        fluence_100 += 1
    if fluence >= 500.0:
        fluence_500 += 1

    width_t = burst.end_acf_range.values.max() - burst.begin_acf_range.values.min()
    ontime = np.sum(burst.end_acf_range.values -  burst.begin_acf_range.values) # comes in ms
    #spectral_lum = speclum(fluence, ontime) / 1e32
    spectral_lum = speclum(fluence, width_t) / 1e32
    scint_bw = burst.scint_bw_MHz.values.astype(float)
    scint_bw_err = burst.scint_bw_err_MHz.values.astype(float)
    mask = (scint_bw == None) | (scint_bw_err == None)
    scint_bw = np.ma.masked_array(scint_bw, mask=mask)
    scint_bw_err = np.ma.masked_array(scint_bw_err, mask=mask)
    weights = 1./np.power(scint_bw_err, 2)
    scint_bw = np.average(scint_bw,
                          weights=weights)
                          #weights=scint_bw_err / scint_bw)
    # scale things to a central freq of 1GHz
    scint_bw *= (1000./cent_f)**4
    scint_bw_err = np.sqrt(1./weights.sum())
    #scint_bw_err = np.sqrt(np.sum(scint_bw_err**2))
    if not scint_bw <= 0.0999:
        scintbws.append(scint_bw)
    # we disregard apparently negative or close-to-zero results
    # for scint_bw (huge errors anyway)
    if (scint_bw <= 0.0999):
        scint_entry = '--$\mathrm{^{h}}$'
    else:
        scint_entry = f'${scint_bw:.1f}\pm{scint_bw_err:.1f}$'
    if (dish == 'St') or (ID == 'B04-o8'):
        row = template1.format(ID, dish, toa, peakSNR, fluence, fluence*0.2, ncomp,
                               width_t, spectral_lum, spectral_lum*0.2,
                               signal_bw, scint_entry)
    else:
        row = template0.format(ID, dish, toa, peakSNR, fluence, fluence*0.2, ncomp,
                               width_t, spectral_lum, spectral_lum*0.2,
                               signal_bw, scint_entry)
    table = table + row

with open('./table.tex', 'w') as f:
    f.write(header)
    f.write(heading)
    f.write('\\hline \n')
    f.write(table)
    f.write(footer)

scintbws = np.array(scintbws)
std = np.std(scintbws)
median = np.median(scintbws)
# from https://stats.stackexchange.com/questions/59838/standard-error-of-the-median
# and https://math.stackexchange.com/questions/2666305/standard-error-of-median#2666367
median_std = 1.253 * std/np.sqrt(len(scintbws))

print(fr'median scintillation bw at 1 GHz: {median} $\pm$ {median_std}')
print(f'We have {fluence_100} bursts above 100 Jyms and {fluence_500} bursts above 500 Jyms.')
