
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from gaiaunlimited.utils import get_healpix_centers
from gaiaunlimited.selectionfunctions import DR3SelectionFunctionTCG, SubsampleSelectionFunction
import healpy as hp


#first, validity of sssf estimate - number of sources in each bin

sssf_df = pd.read_csv(
    '/home/hopkinsl/.gaiaunlimited/hp4_docsRange_rv_mh_fixed.csv',
    comment="#",
    )

print(np.count_nonzero(sssf_df.n<20)/len(sssf_df))

fig,ax = plt.subplots()
ax.hist(sssf_df.n, bins=np.linspace(-0.1,20))
# many have zero - but these could be etremes of colour/mag
# it's fine to have low n if it's because there is actually no stars there
# problem is if I have 'checkerboard' caused by bins being too small, or n=k=1

print(np.count_nonzero((sssf_df.n<20)&(sssf_df.n>0))/len(sssf_df))

fig,ax = plt.subplots()
ax.hist(sssf_df.n, bins=np.linspace(0.5,20))

subsampSF = SubsampleSelectionFunction('', 'hp4_docsRange_rv_mh_fixed',
                           {'healpix': 4,'phot_g_mean_mag': (3,20,0.2), 'g_rp': (-2.5,5.1,0.4)})
#noneed for query, as just loads file

G, col = 17, 0.0
hpc4 = get_healpix_centers(4)
g4 = np.ones_like(hpc4)*G
c = np.ones_like(hpc4)*col

SF, var = subsampSF.query(hpc4, return_variance=True, phot_g_mean_mag_ = g4, g_rp_ = c)
fig = plt.figure(figsize = (20,10))
hp.mollview(SF,fig=fig, hold = True,min = 0,max = 1,title = 'G = {}, G_RP = {}'.format(G,col),coord='CG')



# one and a halfth, average bins over galactic pole to improve 


# hpc = get_healpix_centers(4).galactic

# capI = np.where(np.abs(hpc.b)>30*u.deg)[0]
# plnI = np.where(np.abs(hpc.b)<=30*u.deg)[0]

# print(np.count_nonzero(sssf.iloc[sssf.healpix_].n<20)/len(sssf))

# fig,ax = plt.subplots()
# ax.hist(sssf.iloc[np.logical_not(galCapMask)].n, bins=np.linspace(-0.1,20))

# second, average TCG sf over hp4 pixels
G = 17
survSF = DR3SelectionFunctionTCG()
hpc7 = get_healpix_centers(7)
g7 = np.ones_like(hpc7)*G
SF = survSF.query(hpc7, g7)
fig = plt.figure(figsize = (20,10))
hp.mollview(SF,fig=fig, hold = True,min = 0,max = 1,title = 'G = {}'.format(G),coord='CG')
#dimmest MHRV star is G=16, and urvey is totally complete at G=17, so don't need
#to consider survSF

# from gaiaunlimited.selectionfunctions import DR3SelectionFunctionTCG
# 
# import healpy as hp
# hp.mollview(completeness,coord=['Celestial','Galactic'])
# coords_of_centers = get_healpix_centers(7)
# import numpy as np
# gmag = np.ones_like(coords_of_centers) * 21.
# print(f'Computing the completeness for {len(coords_of_centers)} points.')
# completeness = mapHpx7.query(coords_of_centers,gmag)




# mapHpx7 = DR3SelectionFunctionTCG()













