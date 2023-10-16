import numpy as np

snr = 15
npol = 2
width = 1e-3  # assume a canonocal 1ms-width

order = ['wb_P', 'wb_L', 'st_1', 'st_2', 'O8_1G', 'O8_2G', 'O8_3G', 'Tr_L-2G', 'Tr_L-1G', 'Tr_C-2G', 'Tr_C-1G']
sefds = np.array([2100, 420,  1100, 385,   310,   310,   310,   350,     350,      220,    220])
bw =   np.array([60,   100,   90,   90,   100,   200,    350,  200,      100,      240,   120]) * 1e6

fluence_limits = snr * sefds * np.sqrt(width / (npol * bw)) * 1e3 # in Jy~ms

for band, limit, sefd in zip(order, fluence_limits, sefds):
    print(f'Band: {band}; limit: {limit:0.2f} Jy ms (assumed SEFD: {sefd})')
