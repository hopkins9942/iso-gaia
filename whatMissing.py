import numpy as np
import matplotlib.pyplot as plt
import astropy.io.votable as vt
import healpy as hp
from gaiaunlimited import SubsampleSelectionFunction
from gaiaunlimited.utils import get_healpix_centers

def loadStars(maxD):
    file = vt.parse('/home/hopkinsl/GAIA/data/dataCutProperly-result.vot')
    df = file.get_first_table().to_table(use_names_over_ids=True).to_pandas()
    return df[df.plx>1/maxD]

stars = loadStars(maxD=100)

fig,ax = plt.subplots()
ax.scatter(stars.c, stars.g, s=0.1)

#checking what SSSF misses using SF

SSSF = SubsampleSelectionFunction('', 'avrg_docsRange_rv_mh_age_fixed',
                           {'phot_g_mean_mag': (3,20,0.2), 'g_rp': (-2.5,5.1,0.4)})

sfdf = SSSF.ds.to_dataframe()

misslim=10
print(len(np.where((sfdf.k==0)&(sfdf.n>misslim))[0])/len(sfdf)) # fraction with n but no k

missed = sfdf.iloc[np.where((sfdf.k==0)&(sfdf.n>misslim))[0]]#bins with missed stars

missed_ipix = np.array(missed.index.get_level_values(0))
missed_g = np.array(missed.index.get_level_values(1)) # centre corrds of bins where missed stars
missed_c = np.array(missed.index.get_level_values(2))

fig,ax = plt.subplots()
ax.hist2d(missed_c, missed_g, bins=(np.linspace(missed_c.min()-0.2,
                                    missed_c.max()+0.2,
                                    int((missed_c.max()+0.2-(missed_c.min()-0.2))/0.4)+1),
                                    np.linspace(missed_g.min()-0.1,
                                    missed_g.max()+0.1,
                                    int((missed_g.max()+0.1-(missed_g.min()-0.1))/0.2)+1)))
#missed stars are almost all in color 0-2 and G>13
# spreading of colour at very dim limit could be error
# good - changing mass definitly increases G so as long as effselfunc is never 0 im good




G, col = 17, 0.0
hpc4 = get_healpix_centers(4)
g4 = np.ones_like(hpc4)*G
c = np.ones_like(hpc4)*col

SF, var = SSSF.query(hpc4, return_variance=True, phot_g_mean_mag_ = g4, g_rp_ = c)
fig = plt.figure(figsize = (20,10))
hp.mollview(SF,fig=fig, hold = True,min = 0,max = 1,title = 'G = {}, G_RP = {}'.format(G,col),coord='CG')




# weight = np.zeros()
# for index, vals in stars.iterrows():
#     #get effsel funct
#     weight = 1/effselfunc(vals)
#     isos = (10**feh)*weight
    
#Then bin how I like and plot in histograms

# copositio distribution

# U-V in three composition bins 
# W in comp bins

# speed disrbution to disprove possible ISOs


