import pickle
from math import isclose

import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
import matplotlib.pyplot as plt

import composition as cmp

def arr(gridParams):
    start, stop, step = gridParams
    arr = np.arange(round((stop-start)/step)+1)*step+start
    assert isclose(arr[-1],stop) # will highlight both bugs and when stop-start is not multiple of diff
    return arr

equator = coord.SkyCoord(ra=np.linspace(0,359,360)*u.deg, dec=np.zeros(360)*u.deg)
ecliptic = coord.SkyCoord(lon=np.linspace(0,359,360)*u.deg, lat=np.zeros(360)*u.deg, frame='barycentricmeanecliptic')
galEarthPlane = coord.SkyCoord(l=np.linspace(0,359,360)*u.deg, b=np.zeros(360)*u.deg, frame='galactic')
# galacticPlane = SkyCoord(l=np.linspace(0,359,360)*u.deg, b=np.zeros(360)*u.deg, frame='galactocentric')


#Based on 2k2 solution in Bailer-Jones et al. (2018)
oumuamua = coord.SkyCoord(ra = 279.4752*u.deg,
                             dec = 33.8595*u.deg,
                             distance = 1*u.pc,
                             pm_ra_cosdec = 0*u.mas/u.yr,
                             pm_dec = 0*u.mas/u.yr,
                             radial_velocity = -26.420410*u.km/u.s
                             )

# Based on Bailer-Jones et al. (2020)
borisov = coord.SkyCoord(ra = 32.79720*u.deg,
                             dec = 59.44014*u.deg,
                             distance = 1*u.pc,
                             pm_ra_cosdec = 0*u.mas/u.yr,
                             pm_dec = 0*u.mas/u.yr,
                             radial_velocity = -32.286894*u.km/u.s
                             )

with open('../data/data.pickle', 'rb') as f:
    data = pickle.load(f)

print(data.colnames)
#data does have units, print full col
astrometry = coord.SkyCoord(ra=data['ra'],
                             dec=data['dec'],
                             distance=(data['plx']).to(u.pc, equivalencies=u.parallax()),
                             pm_ra_cosdec=data['pmra'],
                             pm_dec=data['pmdec'],
                             radial_velocity=data['rv']
                             )


# vels = astrometry.velocity

# absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)

MH = data['mh']
fH2O = cmp.comp(MH)

# MH distribution
fig,ax = plt.subplots()
ax.hist(MH, bins=50)
ax.set_xlabel('[M/H]')
fig.savefig('../plots/MH.png', dpi=300)

# fHHO distribution
fig,ax = plt.subplots()
ax.hist(fH2O, bins=50, range=(cmp.fH2OLow, cmp.fH2OHigh))
ax.set_xlabel('fH2O')
ax.set_title('DOESNT INCLUDE MASS/FEH/SLOPE WEIGHTING')
fig.savefig('../plots/fH2O.png', dpi=300)

def angWrap(l):
    return np.where(l>180*u.deg, l-360*u.deg, l)

def spherical(vels):
    #there may be way to use differential object itself to do this
    absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
    theta = np.arcsin(vels.d_z/absv)
    phi = np.arctan2(vels.d_y, vels.d_x)
    return (absv, theta, phi)
    

def velGraph(astrometry, frame='galactic', proj='mollweide', label=''):
    """
    choose frame from ['icrs', 'barycentricmeanecliptic', 'galactic', 'galactocentric', 'galacticlsr', 'lsr']
    for projections see https://matplotlib.org/stable/api/figure_api.html
    mollweide or rectilinear
    """
    
    
    vels = astrometry.transform_to(frame).velocity
    # absv = np.sqrt(vels.d_x**2 + vels.d_y**2 + vels.d_z**2)
    # theta = np.arcsin(vels.d_z/absv)
    # phi = np.arctan2(vels.d_y, vels.d_x)
    absv, theta, phi, = spherical(vels)
    
    
    # projections
    fig, ax = plt.subplots()
    ax.scatter(vels.d_x, vels.d_y, s=0.01, alpha=0.5, color='black')
    ax.scatter(oumuamua.transform_to(frame).velocity.d_x, 
               oumuamua.transform_to(frame).velocity.d_y, s=1, marker='*', color='hotpink')
    ax.scatter(borisov.transform_to(frame).velocity.d_x, 
               borisov.transform_to(frame).velocity.d_y, s=1, marker='*', color='lime')
    ax.set_xlim((-100,100))
    ax.set_ylim((-100,100))
    ax.set_xlabel(r'vx') 
    ax.set_ylabel(r'vy')
    ax.set_aspect('equal')
    ax.set_title(frame+label)
    fig.savefig(f'../plots/vxvy_{frame}_{label}.png', dpi=300)
    
    fig, ax = plt.subplots()
    ax.scatter(vels.d_y, vels.d_z, s=0.01, alpha=0.5, color='black')
    ax.scatter(oumuamua.transform_to(frame).velocity.d_y, 
               oumuamua.transform_to(frame).velocity.d_z, s=1, marker='*', color='hotpink')
    ax.scatter(borisov.transform_to(frame).velocity.d_y, 
               borisov.transform_to(frame).velocity.d_z, s=1, marker='*', color='lime')
    ax.set_xlim((-100,100))
    ax.set_ylim((-30,30))
    ax.set_xlabel(r'vy') 
    ax.set_ylabel(r'vz')
    ax.set_title(frame+label)
    fig.savefig(f'../plots/vyvz_{frame}_{label}.png', dpi=300)
    
    fig, ax = plt.subplots()
    ax.scatter(vels.d_x, vels.d_z, s=0.01, alpha=0.5, color='black')
    ax.scatter(oumuamua.transform_to(frame).velocity.d_x, 
               oumuamua.transform_to(frame).velocity.d_z, s=1, marker='*', color='hotpink')
    ax.scatter(borisov.transform_to(frame).velocity.d_x, 
               borisov.transform_to(frame).velocity.d_z, s=1, marker='*', color='lime')
    ax.set_xlim((-100,100))
    ax.set_ylim((-30,30))
    ax.set_xlabel(r'vx') 
    ax.set_ylabel(r'vz')
    ax.set_title(frame+label)
    fig.savefig(f'../plots/vxvz_{frame}_{label}.png', dpi=300)
    
    
    # on-sky
    fig = plt.figure()
    ax = fig.add_subplot(projection=proj)
    ax.scatter(phi.to(u.deg), theta.to(u.deg), s=0.01, alpha=0.5, color='black')
    # lineS = 0.5
    # if (frame=='galactic')or(frame=='galacticlsr'):
    #     ax.scatter(angWrap(equator.galactic.l).to(u.rad), equator.galactic.b.to(u.rad), s=lineS, color='hotpink')
    #     ax.scatter(angWrap(ecliptic.galactic.l).to(u.rad), ecliptic.galactic.b.to(u.rad), s=lineS, color='lime')
    # if (frame=='icrs')or(frame=='lsr'):
    #     ax.scatter(angWrap(ecliptic.icrs.ra).to(u.rad), ecliptic.icrs.dec.to(u.rad), s=lineS, color='lime')
    #     ax.scatter(angWrap(galEarthPlane.icrs.ra).to(u.rad), galEarthPlane.icrs.dec.to(u.rad), s=lineS, color='cyan')
    # if frame=='barycentricmeanecliptic':
    #     ax.scatter(angWrap(equator.barycentricmeanecliptic.lon).to(u.rad), equator.barycentricmeanecliptic.lat.to(u.rad), s=lineS, color='hotpink')
    #     ax.scatter(angWrap(galEarthPlane.barycentricmeanecliptic.lon).to(u.rad), galEarthPlane.barycentricmeanecliptic.lat.to(u.rad), s=lineS, color='cyan')
    ax.scatter(spherical(oumuamua.transform_to(frame).velocity)[2].to(u.deg), 
               spherical(oumuamua.transform_to(frame).velocity)[1].to(u.deg), s=1, marker='*', color='hotpink')
    ax.scatter(spherical(borisov.transform_to(frame).velocity)[2].to(u.deg), 
               spherical(borisov.transform_to(frame).velocity)[1].to(u.deg), s=1, marker='*', color='lime')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\theta$')
    ax.set_xlim(180, -180)
    ax.set_ylim(-90, 90)
    # can't reverse axes of mollweide, try cartopy
    ax.set_title(frame+label)
    fig.savefig(f'../plots/onsky_{frame}_{label}.png', dpi=300)
    
    
    # absv
    fig,ax = plt.subplots()
    ax.hist(absv, bins=50)
    ax.set_xlabel(r'speed/ km/s')
    # ax.set_xlim((0,150))
    ax.set_title(frame+label)
    fig.savefig(f'../plots/speed_{frame}_{label}.png', dpi=300)
    

if __name__=='__main__':
    # binEdges = np.array([cmp.FeHLow+0.0001, cmp.compInv(0.3), cmp.FeHHigh-0.0001])
    # labels = ['low', 'midlow', 'midhigh', 'high'][::-1]
    
    binEdges = np.array([0])
    labels = ['subsolar', 'supersolar']
    
    bindex = np.digitize(MH, binEdges) # index of mh bin
    frameNames = ['galactic', 'galacticlsr', 'galactocentric']
    for i in range(len(binEdges)+1):
        for f in frameNames:
            velGraph(astrometry[bindex==i], frame=f, label=labels[i], proj='rectilinear')





