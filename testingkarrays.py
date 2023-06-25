import numpy as np
from astropy.io import votable as vt
import pandas as pd
import healpy as hp



file_apinside = vt.parse('/home/hopkinsl/GAIA/data/just_k_outsidejoin_apinside-result.vot')
df_apinside = file_apinside.get_first_table().to_table(use_names_over_ids=True).to_pandas()

file_outsidejoin = vt.parse('/home/hopkinsl/GAIA/data/just_k_outsidejoin-result.vot')
df_outsidejoin = file_outsidejoin.get_first_table().to_table(use_names_over_ids=True).to_pandas()

file_subjoin = vt.parse('/home/hopkinsl/GAIA/data/just_k_subjoin-result.vot')
df_subjoin = file_subjoin.get_first_table().to_table(use_names_over_ids=True).to_pandas()

# al same length
print(len(df_apinside),
      len(df_outsidejoin),
      len(df_subjoin),
      )

#diffs

print((df_apinside-df_outsidejoin).max())
print((df_apinside-df_subjoin).max())

print((df_apinside-df_outsidejoin).min())
print((df_apinside-df_subjoin).min())

"""All checks out! Use apinside, it's faster o query on gaia"""