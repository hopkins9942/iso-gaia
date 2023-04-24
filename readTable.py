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
with open('../data/data.pickle', 'rb') as f:
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

print(astrometryICRS.velocity)
#astrometryLSR = astrometryICRS.transform_to(coord.LSR)
#print(astrometryLSR)

# print(coord.Galactocentric.galcen_v_sun)

vels = astrometryICRS.velocity

absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)

fig,ax = plt.subplots()
ax.hist(absv[lowMH], bins=50, alpha=0.5, label='lowMH')
ax.hist(absv[highMH], bins=50, alpha=0.5, label='highMH')
ax.set_xlabel(r'speed/ km/s')
ax.legend()
fig.savefig('../speed.png', dpi=300)

theta = np.arcsin(vels.d_z/absv)#.to(u.deg)
phi = np.arctan2(vels.d_y, vels.d_x)#.to(u.deg)

fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(phi[lowMH], theta[lowMH], s=0.1, alpha=0.5, label='lowMH')
ax.scatter(phi[highMH], theta[highMH], s=0.1, alpha=0.5, label='highMH')
# ax.legend()
ax.set_xlabel(r'$l$')
ax.set_ylabel(r'$b$')
fig.savefig('../angle.png', dpi=300)


