# Code to figure out at what maximum redshift an FRB would still be
# seen by FAST if we detect it with a certain fluence and a certain distance
# at one of our smaller dishes.

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import Planck15

fluence = 1593 * 0.8 # 1.2 #0.8              # Jy ms
z_source = 0.098            # for R67
a = 0.0                     # because we give equal weights to all nu when computing the fluence
BW = 90e6                  # effective BW at Stockert
detection_threshold = 53e-3 # Jy ms, from FAST paper

clancy = False

if clancy:
    z_source = 0.522
    a = -1.5
    BW = 1e9
    detection_threshold = 1

lum_dist = (Planck15.luminosity_distance(z_source)).to(u.cm).value

E = fluence * 1e-3 * 4 * np.pi * lum_dist * lum_dist * BW * 1e-23 / (1+z_source)**(2+a)

print(f'Energy of the burst: {np.log10(E):.2f}')

if clancy:
    E = 10**(41.7)

# redshift range to probe
z = np.linspace(0.1, 20.0, 1000)

# gotta take care of spectral index here
a = -1.5                   # from Clancy James paper https://ui.adsabs.harvard.edu/abs/2022MNRAS.510L..18J/abstract

# we try and find where the luminosity distance that we get from the energy for a certain redshift
# equals that of what Planck15 gets for that redshift.
ratio = np.sqrt((E * (1+z)**(2+a))/ (detection_threshold * 1e-3 * 4 * np.pi * BW * 1e-23)) / (Planck15.luminosity_distance(z)).to(u.cm).value

bestz = z[np.argmin(abs(ratio-1))]
print(f'Best z is {bestz}.')
# plot it up
plt.plot(z, ratio);
plt.hlines(1, z.min(), z.max())
plt.xlabel('z')
plt.ylabel('D_L(E, z) / D_L(z) per Planck15')
plt.ylim([0, 2])
plt.show()
