#!/usr/bin/env python3
import numpy as np

print('Volumentric rates based on energy distro\n')
# Shin et al. 2023 volumetric rates above pivot energy E_p = 1e39erg from CHIME/FRB cat
E_p = 1e39          # erg, pivot energy at 600 MHz
chime_rate = 7.3e4  # per GPc per year at 600 MHz; +8.8-3.8 above E_p
bw = 300e6          # fiducial bandwidth of bursts, MHz
alpha = -1.5        # spectral index
gamma = -0.3        # slope of cummulative energy distro (in Shin et al. it's differential
nu_chime = 600e6    # CHIME frequency MHz
nu_precise = 1300e6 # our central frequency MHz
n_bursts = 42       # number of bursts detected above 3e30 erg/Hz;
E_nu_max = 1.3e32   # brightest end; erg/Hz
n_brightest = 3     # number of bursts above 1.3e32erg/Hz (i.e. > 500 Jyms)
hrs_observed = 2281 # total number of non-overlapping hours on source
distance = 0.453    # Gpc; luminosity distance to R67

# spectral index makes the same bursts fainter but we keep the pivot energy constant at 1e39; so a burst that has >1e39 at 1300 MHz must be
# accordingly brighter at 600 MHz. So what's the rate of those bursts at 600 MHz?
E_600 = E_p * (600. / 1300.)**(alpha)           # energy of burst at 600 MHz that has an energy of E_p at 1300 MHz
rate_600 = chime_rate * (E_600 / E_p)**(gamma)  # volumetric rate

print(f'Volumetric rate at 1300 MHz is {round(rate_600,1)}/Gpc/yr.')

E_nu_p_600 = E_p / bw  # convert to spectral energy; erg/Hz

volume = 4./3. * np.pi * distance**3
fractional_year = hrs_observed / 8760.

# so we observed n_bursts within volume during fractional_year; this corresponds to
observed_fraction = n_bursts * 1./volume * 1./fractional_year

print(f'Given the distance, time on source and volume covered, this corresponds to {round(observed_fraction,2)}/Gpc/yr')
print(f'or {round(observed_fraction/rate_600*100.,2)}%')

# so we observed n_brightest within volume during fractional_year; this corresponds to
observed_fraction = n_brightest * 1./volume * 1./fractional_year

# and we expect a volumetric rate of
E_nu_600_max = E_nu_max * (600. / 1300.)**(alpha) # energy of bursts at 600 MHz that have spectral energy of 1.3e32 erg/Hz at 1300 MHz
rate_max = chime_rate * (E_nu_600_max / (E_p/bw))**(gamma)
print(f'Considering the 3 brightest bursts only this corresponds to {round(observed_fraction,2)}/Gpc/yr')
print(f'or {round(observed_fraction/rate_max*100.,2)}%')


print('\n Skyrates based on fluence distro\n')
# N(>F) = 20 * (F / 50 Jyms)**(-1.5); from Lu et al. 2019
total_exposure = 5.1e5  # deg^2*hr; ASKAP survey, Shannon et al. 2018, Nature
skyrate_100 = 5e3  # /sky/yr above 100 Jyms
skyrate_500 = skyrate_100 / (1./5.)**(-1.5)
nbursts_100 = 17  # bursts with fluence >= 100 Jyms (from pandas2latex)
nbursts_500 = 3   # burst with fluence >= 500 Jyms (from pandas2latex)

print(f'We have {round(nbursts_100 / fractional_year,1)} bursts/yr above 100 Jyms and {round(nbursts_500 / fractional_year,1)} bursts/yr above 500 Jyms')
print(f'That corresponds to {round(nbursts_100 / fractional_year / skyrate_100 * 100,1)} and {round(nbursts_500 / fractional_year / skyrate_500 * 100,1)} %, respectivley.')
