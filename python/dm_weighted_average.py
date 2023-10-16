import numpy as np

dms = {
    'p67113_1':[410.833, 0.335],
    'p67113_2':[410.621, 0.580],
    'p67120_2':[410.735, 0.558],
    'p67125_2':[411.053, 0.265],
    'p67138_2':[410.408, 0.426],
    'pre028':[411.334, 0.341],
    'r67022':[411.886, 0.277],
    'r67028':[410.935, 0.253],
    'r67079_2':[411.189, 0.350],
    'r67081_1':[410.774, 0.286]
}

dm = []
dm_errs = []
for key, val in dms.items():
    dm.append(val[0])
    dm_errs.append(val[1])

dm = np.array(dm)
dm_errs = np.array(dm_errs)
weights = 1./np.power(dm_errs, 2)

av = np.average(dm, weights=weights)
err = np.sqrt(1./weights.sum())

print(f'DM = {av} +/- {err}')
