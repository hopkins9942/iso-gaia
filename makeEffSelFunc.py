import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.votable as vt
import healpy as hp
from gaiaunlimited import SubsampleSelectionFunction

import myIsochrones


"""Bin in D, age and FeH. calculate ESF in each assume averaged over sky,
 so can for each star use value of bin it's in 
"""

sssf = pd.read_csv('/home/hopkinsl/GAIA/data/sssf.csv', comment="#")

print(sssf)
ss_shape = (85,19)
sssf_array = np.zeros(ss_shape)
for index, vals in sssf.iterrows():
    g_i, c_i, n, k, p = vals
    if n>=5: 
        sssf_array[int(g_i), int(c_i)] = p #finally in a form I'm used to using
print(sssf_array)

fig,ax = plt.subplots()
ax.imshow(sssf_array.T, origin='lower')

#this if fine as long as isochrones don't sproduce stars that don't exist in universe,
#meaning where n=k=0

G_edges = np.linspace(3,20,85+1)
c_edges = np.linspace(-2.5, 5.1, 19+1)


e_shape = (10, 15, 14)
D_edges = np.linspace(0, 200, 10+1)
mh_edges = np.linspace(-1, 0.5, 15+1)
age_edges = np.linspace(0, 14, 14+1)

isogrid = myIsochrones.loadGrid()
MH_logAge, indx = myIsochrones.extractIsochrones(isogrid)
indx = np.append(indx, len(isogrid))

Mmax = np.array([isogrid['Mini'][indx[i+1]-2]] for i in range(len(indx)-1))
# indx[i+1] is start of i+1th isochrone. indx[i+1]-2 is second last point of ith, ie highest mass non-wd
D_mid = (D_edges[:-1]+D_edges[1:])/2
mus = 5*(np.log(D_mid)-1)
esf_array = np.zeros(e_shape)
for i in range(len(MH_logAge)):
    isochrone = isogrid[indx[i]:indx[i+1]-1]#skips WD at end
    weights = myIsochrones.calcWeights(isochrone)
    mh_i, age_i = np.unravel_index(i, e_shape[1:])
    
    # print(MH_logAge[i][0],MH_logAge[i][1])
    # print(mh_i, age_i)
    c_array = isochrone['Gmag'] - isochrone['G_RPmag']
    # print(c_array)
    c_indices = np.clip(np.digitize(c_array, c_edges)-1, 0, len(c_edges)-2)
    # print(c_indices)
   
    
    for D_i, mu in enumerate(mus):
        G_array = isochrone['Gmag'] + mu
        G_indices = np.clip(np.digitize(G_array, G_edges)-1, 0, len(G_edges)-2)
        for isopoint in range(len(isochrone)):
            esf_array[D_i, mh_i, age_i]+= sssf_array[G_indices[isopoint], c_indices[isopoint]]*weights[isopoint]
            
assert np.count_nonzero(esf_array==0) == 0 # yeah boiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

fig,ax = plt.subplots()
ax.imshow(esf_array[5,:,:].T, origin='lower')
#hmm

# print(MH_logAge)
# print(indx)




