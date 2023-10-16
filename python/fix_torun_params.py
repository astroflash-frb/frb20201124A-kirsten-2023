import pandas as pd

sefd = 350.    # the new value for the SEFD
station = 'tr' # the station to fix
db = '../dbs/burst_info.pickle' # pandas dataframe that carries the data

cols2fix = ['fluence_jyms', 'peak_flux_jy', 'energy_iso_erg_per_hz',
            'spectral_lum_erg_per_s_per_hz']

df = pd.read_pickle(db)

sefd_old = float(df[(df.dish == station)].sefd_jy.unique()[0])

ratio = sefd / sefd_old

# *should* be just a linear scaling, famous last words.
for col in cols2fix:
    df.loc[(df.dish == station), col] *= ratio

# also fix the SEFD column
df.loc[(df.dish == station), 'sefd_jy'] = int(sefd)

df.to_pickle(db)
