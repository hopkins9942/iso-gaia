import numpy as np
from astropy.io import votable as vt
import pandas as pd
import healpy as hp

"""Assumes n and k sel funct arrays are made externally and cmbines, 
 output is saved file that can be picked p and used by
gaiaunlimited subsampleselectionfunction with missing bins added"""

print('start')
k_file = vt.parse('/home/hopkinsl/GAIA/data/just_k_rv-mh-age-result.vot')
k_df = k_file.get_first_table().to_table(use_names_over_ids=True).to_pandas()

print('k loaded')
n_file = vt.parse('/home/hopkinsl/GAIA/data/just_n-result.vot')
n_df = n_file.get_first_table().to_table(use_names_over_ids=True).to_pandas()

print('n loaded')
print(k_df)
print(n_df)

# print('looking for oddities')
# for tempdf in [k_df,n_df]:
#     print(tempdf.healpix_.min())
#     print(tempdf.healpix_.max())
#     print(tempdf.phot_g_mean_mag_.min())
#     print(tempdf.phot_g_mean_mag_.max())
#     print(tempdf.g_rp_.min())
#     print(tempdf.g_rp_.max())

# print('removing oddities')
# k_df = k_df[k_df.g_rp_!=18.0]
# n_df = n_df[n_df.g_rp_!=18.0]

# assert len(k_df)==255200
# assert len(n_df)==1047689 - 6008

shape = (hp.order2npix(4), 85, 19) 
print(shape)
n_col = np.zeros(np.prod(shape))
k_col = np.zeros(np.prod(shape))


print('starting k scan')
for i, row in k_df.iterrows():
    hpi, gi, ci, k = row
    index = np.ravel_multi_index((int(hpi),int(gi),int(ci)), shape)
    k_col[index] = k
    
print('starting n scan')
for i, row in n_df.iterrows():
    hpi, gi, ci, n = row
    index = np.ravel_multi_index((int(hpi),int(gi),int(ci)), shape)
    n_col[index] = n

print('finished scan')
hpi_col, gi_col, ci_col = np.unravel_index(np.arange(np.prod(shape)), shape)

datadict = {'healpix_':hpi_col,
            'phot_g_mean_mag_':gi_col,
            'g_rp_':ci_col,
            'n':n_col,
            'k':k_col}

print('making nk_df')
nk_df = pd.DataFrame(datadict)

hplevel_and_binning = {'healpix': 4,'phot_g_mean_mag': (3,20,0.2), 'g_rp': (-2.5,5.1,0.4)}

print('writing')
with open('/home/hopkinsl/.gaiaunlimited/hp4_docsRange_rv_mh_age_fixed.csv', "w") as f:
    f.write(f"#{hplevel_and_binning}\n")
    nk_df.to_csv(f, index = False)

print('done!')
# sssf_array = (k_array+1)/(n_array+2)
# var_array = sssf_array*(n_array-k_array+1)/((n_array+2)*(n_array+3))



