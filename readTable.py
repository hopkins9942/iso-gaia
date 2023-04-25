import os
import time
import numpy as np
import pickle
from math import isclose

import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from astropy.table import Table

import astropy.units as u
from astropy.coordinates import SkyCoord


def arr(gridParams):
    start, stop, step = gridParams
    arr = np.arange(round((stop-start)/step)+1)*step+start
    assert isclose(arr[-1],stop) # will highlight both bugs and when stop-start is not multiple of diff
    return arr

equator = SkyCoord(ra=np.linspace(0,359,360)*u.deg, dec=np.zeros(360)*u.deg)
ecliptic = SkyCoord(lon=np.linspace(0,359,360)*u.deg, lat=np.zeros(360)*u.deg, frame='barycentricmeanecliptic')
galEarthPlane = SkyCoord(l=np.linspace(0,359,360)*u.deg, b=np.zeros(360)*u.deg, frame='galactic')
# galacticPlane = SkyCoord(l=np.linspace(0,359,360)*u.deg, b=np.zeros(360)*u.deg, frame='galactocentric')

with open('../data/data.pickle', 'rb') as f:
    table = pickle.load(f)

print(table.colnames)

astrometryICRS = coord.SkyCoord(ra=table['ra'],
                             dec=table['dec'],
                             distance=(table['plx']).to(u.pc, equivalencies=u.parallax()),
                             pm_ra_cosdec=table['pmra'],
                             pm_dec=table['pmdec'],
                             radial_velocity=table['rv']
                             )
#astrometryLSR = astrometryICRS.transform_to(coord.LSR)
#print(astrometryLSR)

# print(coord.Galactocentric.galcen_v_sun)

vels = astrometryICRS.velocity

absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)

binEdges = arr((-1,0.5,0.5))
# labels = ['vlow', 'low', 'mid', 'high', 'vhigh']
bindex = np.digitize(table['mh'], binEdges) # index of mh bin

# MH distribution
fig,ax = plt.subplots()
ax.hist(table['mh'], bins=50)
ax.set_xlabel('[M/H]')
fig.savefig('../plots/MH.png', dpi=300)

# speed
fig,ax = plt.subplots()
for i in range(len(binEdges)+1):
    ax.hist(absv[bindex==i], bins=50, alpha=0.5, label=f'{i}')
ax.set_xlabel(r'speed/ km/s')
ax.set_xlim((0,150))
ax.legend()
fig.savefig('../plots/speed.png', dpi=300)


theta = np.arcsin(vels.d_z/absv)#.to(u.deg)
phi = np.arctan2(vels.d_y, vels.d_x)#.to(u.deg)

fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
for i in range(len(binEdges)+1):
    ax.scatter(phi[bindex==i], theta[bindex==i], s=0.1, alpha=0.5, label=f'{i}')
ax.legend()
# ax.set_xticklabels([])
# ax.set_yticklabels([])
ax.set_xlabel(r'ra') #are these galactic coordinates? ICRS?
ax.set_ylabel(r'dec')
fig.savefig('../plots/angle.png', dpi=300)

# trying other coordinate systems
vels = astrometryICRS.galactic.velocity

absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
theta = np.arcsin(vels.d_z/absv)#.to(u.deg)
phi = np.arctan2(vels.d_y, vels.d_x)#.to(u.deg)
fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
for i in range(len(binEdges)+1):
    ax.scatter(phi[bindex==i], theta[bindex==i], s=0.1, alpha=0.5, label=f'{i}')
ax.legend()
# ax.set_xticklabels([])
# ax.set_yticklabels([])
ax.set_xlabel(r'l') #are these galactic coordinates? ICRS?
ax.set_ylabel(r'b')
fig.savefig('../plots/angle.png', dpi=300)


# trying other coordinate systems
vels = astrometryICRS.galacticlsr.velocity

absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
theta = np.arcsin(vels.d_z/absv)#.to(u.deg)
phi = np.arctan2(vels.d_y, vels.d_x)#.to(u.deg)
fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
for i in range(len(binEdges)+1):
    ax.scatter(phi[bindex==i], theta[bindex==i], s=0.1, alpha=0.5, label=f'{i}')
# ax.legend()
# ax.set_xticklabels([])
# ax.set_yticklabels([])
ax.set_xlabel(r'l') #are these galactic coordinates? ICRS?
ax.set_ylabel(r'b')
fig.savefig('../plots/angle.png', dpi=300)

def wrapper(l):
    return np.where(l>180*u.deg, l-360*u.deg, l)

def velPlots(sc, frame, proj='mollweide'):
    """
    choose icrs, galactic, galacticlsr for frame
    for projections see https://matplotlib.org/stable/api/figure_api.html
    mollweide or rectilinear
    """
    
    vels = sc.transform_to(frame).velocity
    absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
    theta = np.arcsin(vels.d_z/absv).to(u.deg)
    phi = np.arctan2(vels.d_y, vels.d_x).to(u.deg)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection=proj)
    ax.scatter(phi, theta, s=0.1, alpha=0.5, color='grey')
    lineS = 0.5
    if (frame=='galactic')or(frame=='galacticlsr'):
        ax.scatter(wrapper(equator.galactic.l), equator.galactic.b, s=lineS, color='hotpink')
        ax.scatter(wrapper(ecliptic.galactic.l), ecliptic.galactic.b, s=lineS, color='lime')
    if (frame=='icrs')or(frame=='lsr'):
        ax.scatter(wrapper(ecliptic.icrs.ra), ecliptic.icrs.dec, s=lineS, color='lime')
        ax.scatter(wrapper(galEarthPlane.icrs.ra), galEarthPlane.icrs.dec, s=lineS, color='cyan')
    if frame=='barycentricmeanecliptic':
        ax.scatter(wrapper(equator.barycentricmeanecliptic.lon), equator.barycentricmeanecliptic.lat, s=lineS, color='hotpink')
        ax.scatter(wrapper(galEarthPlane.barycentricmeanecliptic.lon), galEarthPlane.barycentricmeanecliptic.lat, s=lineS, color='cyan')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\theta$')
    ax.set_title(frame)
    
    fig, ax = plt.subplots()
    ax.scatter(vels.d_x, vels.d_y, s=0.1, alpha=0.5, color='black')
    # ax.set_xlim(())
    ax.set_xlabel(r'vx') 
    ax.set_ylabel(r'vy')
    ax.set_title(frame)
    
    fig, ax = plt.subplots()
    ax.scatter(vels.d_y, vels.d_z, s=0.1, alpha=0.5, color='black')
    # ax.legend()
    ax.set_xlabel(r'vy') 
    ax.set_ylabel(r'vz')
    ax.set_title(frame)
    
    fig, ax = plt.subplots()
    ax.scatter(vels.d_x, vels.d_z, s=0.1, alpha=0.5, color='black')
    # ax.legend()
    ax.set_xlabel(r'vx') 
    ax.set_ylabel(r'vz')
    ax.set_title(frame)

for f in ['icrs', 'barycentricmeanecliptic', 'galactic', 'galactocentric', 'galacticlsr', 'lsr']:
    velPlots(astrometryICRS, f, proj='rectilinear')
    
    
    
    
    
    
    



# def velPlots(sc, frame, proj='mollweide'):
#     """
#     choose icrs, galactic, galacticlsr for frame
#     for projections see https://matplotlib.org/stable/api/figure_api.html
#     mollweide or rectilinear
#     """
    
#     vels = sc.transform_to(frame).velocity
#     absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
#     theta = np.arcsin(vels.d_z/absv)
#     phi = np.arctan2(vels.d_y, vels.d_x)
    
#     fig = plt.figure()
#     ax = fig.add_subplot(projection=proj)
#     for i in range(len(binEdges)+1):
#         ax.scatter(phi[bindex==i], theta[bindex==i], s=0.1, alpha=0.5, label=f'{i}')
#     if 'galact' in frame:
#         ax.plot
#     # ax.legend()
#     # ax.set_xticklabels([])
#     # ax.set_yticklabels([])
#     ax.set_xlabel(r'$\phi$') #are these galactic coordinates? ICRS?
#     ax.set_ylabel(r'$\theta$')
#     ax.set_title(frame)
    
#     fig, ax = plt.subplots()
#     for i in range(len(binEdges)+1):
#         ax.scatter(vels.d_x, vels.d_y, s=0.1, alpha=0.5, label=f'{i}')
#     # ax.legend()
#     # ax.set_xlim(())
#     ax.set_xlabel(r'vx') 
#     ax.set_ylabel(r'vy')
#     ax.set_title(frame)
    
#     fig, ax = plt.subplots()
#     for i in range(len(binEdges)+1):
#         ax.scatter(vels.d_y, vels.d_z, s=0.1, alpha=0.5, label=f'{i}')
#     # ax.legend()
#     ax.set_xlabel(r'vy') 
#     ax.set_ylabel(r'vz')
#     ax.set_title(frame)
    
#     fig, ax = plt.subplots()
#     for i in range(len(binEdges)+1):
#         ax.scatter(vels.d_x, vels.d_z, s=0.1, alpha=0.5, label=f'{i}')
#     # ax.legend()
#     ax.set_xlabel(r'vx') 
#     ax.set_ylabel(r'vz')
#     ax.set_title(frame)
    

    
    
    
    
    
    
    
    


