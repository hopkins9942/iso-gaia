import numpy as np
import matplotlib.pyplot as plt
import astropy.io.votable as vt
import healpy as hp
from gaiaunlimited import SubsampleSelectionFunction


def loadStars(maxD):
    file = vt.parse('/home/hopkinsl/GAIA/data/dataCutProperly-result.vot')
    df = file.get_first_table().to_table(use_names_over_ids=True).to_pandas()
    return df[df.plx>1000/maxD]

# def effSelFunc()


stars = loadStars(maxD=200)
print(len(stars))
print(min(stars.mh))
print(max(stars.mh))
print(min(stars.age))
print(max(stars.age))

#check roughly uniform in sky and dist
fig,ax = plt.subplots()
ax.scatter(stars.ra, stars.dec, s=0.1)

fig,ax = plt.subplots()
ax.hist(1/stars.plx, bins=50)
ax.plot(np.linspace(0,0.1), 3600*(np.linspace(0,0.1)/0.1)**2)
ax.plot(np.linspace(0,0.2), 9000*(np.linspace(0,0.2)/0.2))

#not unifor in D - edge of disk visible at 200pc
# However, still looks fairly D**2 at 100pc - good enough for this
# Acually don't need to assume uniform in D - calculate ESF in M_G, age and mh
# This for each 

#Mking effsel func

# SSSF = SubsampleSelectionFunction('', 'avrg_docsRange_rv_mh_age_fixed',
#                            {'phot_g_mean_mag': (3,20,0.2), 'g_rp': (-2.5,5.1,0.4)})

# sfdf = SSSF.ds.to_dataframe()