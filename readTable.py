import os
import time
import numpy as np
import pickle

import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from astropy.table import Table

import astropy.units as u
from astropy.coordinates import SkyCoord


# table= Table.read('data.dat', format='ascii')
with open('data.pickle', 'rb') as f:
    table = pickle.load(f)


print(table)

lowMH = table['mh']<0
highMH = table['mh']>=0

print(lowMH)
print(highMH)

fig,ax = plt.subplots()
ax.hist(abs(table[lowMH]['rv']), bins=50, alpha=0.5, label='lowMH')
ax.hist(abs(table[highMH]['rv']), bins=50, alpha=0.5, label='highMH')
ax.set_xlabel('Radial Velocity / km/s')
ax.set_xlim((0,150))
ax.legend()

astrometryICRS = coord.SkyCoord(ra=table['ra'],
                             dec=table['dec'],
                             distance=(table['plx']).to(u.pc, equivalencies=u.parallax()),
                             pm_ra_cosdec=table['pmra'],
                             pm_dec=table['pmdec'],
                             radial_velocity=table['rv']
                             )
#units in table?
print(astrometryICRS)

astrometryLSR = astrometryICRS.transform_to(coord.LSR)
print(astrometryLSR)

# print(coord.Galactocentric.galcen_v_sun)

# v = np.sqrt(astrometryGC.v_x**2 + astrometryGC.v_y**2 + astrometryGC.v_z**2)

# fig,ax = plt.subplots()
# ax.hist(v[lowMH], bins=50, alpha=0.5, label='lowMH')
# ax.hist(v[highMH], bins=50, alpha=0.5, label='highMH')
# ax.set_xlabel('Radial Velocity / km/s')
# ax.legend()

# theta = np.arccos(astrometryGC.v_z/v)
# phi = np.arctan2(astrometryGC.v_y, astrometryGC.v_x)

# fig,ax = plt.subplots()
# ax.scatter(phi[lowMH], theta[lowMH], s=0.1, label='lowMH')
# ax.scatter(phi[highMH], theta[highMH], s=0.1,  label='highMH')
# ax.legend()
