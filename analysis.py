import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.votable as vt
import astropy.units as u
import astropy.coordinates as coord
import healpy as hp
from gaiaunlimited import SubsampleSelectionFunction


def loadStars(maxD):
    file = vt.parse('/home/hopkinsl/GAIA/data/dataCutProperly-result.vot')
    df = file.get_first_table().to_table(use_names_over_ids=True).to_pandas()
    return df[df.plx>1000/maxD]


# Defining comp:
FeH_p = np.array([-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4])
fH2O_p = np.array([0.5098, 0.4905, 0.4468, 0.4129, 0.3563, 0.2918, 0.2173, 0.1532, 0.06516])
compPoly = np.polynomial.polynomial.Polynomial.fit(FeH_p, fH2O_p, 3)
FeHLow = FeH_p[0]
FeHHigh = FeH_p[-1]
fH2OLow = compPoly(FeH_p[-1])
fH2OHigh = compPoly(FeH_p[0])
def comp(FeH):
    return np.where(FeHLow<=FeH, np.where(FeH<FeHHigh,
                                          compPoly(FeH),
                                          fH2OLow), fH2OHigh)
    
def compInv(fH2O):
    """inv may not work with array inputs"""
    if np.ndim(fH2O)==0:
        val = fH2O
        if fH2OLow<=val<=fH2OHigh:
            allroots = (compPoly-val).roots()
            myroot = allroots[(FeH_p[0]<=allroots)&(allroots<=FeH_p[-1])]
            assert len(myroot)==1 # checks not multiple roots
            assert np.isreal(myroot[0])
            return np.real(myroot[0])
        else:
            return np.nan
    else:
        returnArray = np.zeros_like(fH2O)
        for i, val in enumerate(fH2O):
            if fH2OLow<=val<fH2OHigh:
                allroots = (compPoly-val).roots()
                myroot = allroots[(FeH_p[0]<=allroots)&(allroots<FeH_p[-1])]
                assert len(myroot)==1 # checks not multiple roots
                assert np.isreal(myroot[0])
                returnArray[i] = np.real(myroot[0])
            else:
                returnArray[i] =  np.nan
        return returnArray
        
def compDeriv(FeH):
    return np.where((FeHLow<=FeH)&(FeH<FeHHigh), compPoly.deriv()(FeH), 0)

for x in [-0.4+0.0001, -0.2, 0, 0.2, 0.4-0.0001]:
    assert np.isclose(compInv(comp(x)), x) #checks inverse works
    

# if exists load else copute 
# stars = loadStars(maxD=200)
with open('/home/hopkinsl/GAIA/data/stars200.pickle', 'rb') as f:
    stars = pickle.load(f)
print(stars.columns)
stars = stars.to_numpy() #so Iknow what I' doing - may switch back when back on internet
# for now use 2 for plx, 4-6 for ra,dec,pmra, pmdec, 7 for rv, 8-9 for mh,age
assert stars.shape==(197243,11)

with open('/home/hopkinsl/GAIA/data/esf.pickle', 'rb') as f:
    esf_array = pickle.load(f)

def ESFquery(D, mh, age):
    D_edges = np.linspace(0, 200, 10+1)
    mh_edges = np.linspace(-1, 0.5, 15+1)
    age_edges = np.linspace(0, 14, 14+1)
    D_i = np.clip(np.digitize(D, D_edges)-1, 0, len(D_edges)-2)
    mh_i = np.clip(np.digitize(mh, mh_edges)-1, 0, len(mh_edges)-2)
    age_i = np.clip(np.digitize(age, age_edges)-1, 0, len(age_edges)-2)
    
    return esf_array[D_i, mh_i, age_i]

# starWeights = np.zeros(197243)
# for i in range(197243):
#     starWeights[i]=1/ESFquery(1000/stars[i,2], stars[i,8], stars[i,9])
    
# with open('/home/hopkinsl/GAIA/data/starWeights200.pickle', 'wb') as f:
#     pickle.dump(starWeights,f)

with open('/home/hopkinsl/GAIA/data/starWeights200.pickle', 'rb') as f:
    starWeights = pickle.load(f) 

# starcoords = coord.SkyCoord(ra=stars.ra.to_numpy()*u.deg,
#                              dec=stars.dec.to_numpy()*u.deg,
#                               distance=(1000/stars.plx.to_numpy())*u.pc,
#                               pm_ra_cosdec=stars.pmra.to_numpy()*u.mas/u.yr,
#                               pm_dec=stars.pmdec.to_numpy()*u.mas/u.yr,
#                               radial_velocity=stars.rv.to_numpy()*u.km/u.s
#                              )#units don't seem to work with pd.Series

# starcoords = coord.SkyCoord(ra=stars[:,3]*u.deg,
#                               dec=stars[:,4]*u.deg,
#                               distance=(1000/stars[:,2])*u.pc,
#                               pm_ra_cosdec=stars[:,5]*u.mas/u.yr,
#                               pm_dec=stars[:,6]*u.mas/u.yr,
#                               radial_velocity=stars[:,7]*u.km/u.s
#                               )
# vels = starcoords.transform_to('galactic').velocity
# print(vels)
# with open('/home/hopkinsl/GAIA/data/starVels200.pickle', 'wb') as f:
#     pickle.dump(vels,f)

with open('/home/hopkinsl/GAIA/data/starVels200.pickle', 'rb') as f:
    vels = pickle.load(f)


# isos = np.zeros((stars.shape[0], 5)) # comp, U,V,W, weight - rubbish but will have to do till Im back on internet
# # i = 0
# for i in range(197243):
#     # if i%1000==0:
#     #     print(i) doesn't work as i is index, which is kept from whoe file so not complete
#     isos[i,0] = comp(stars[i,8])
#     isos[i,1] =  vels[i].d_x/(u.km/u.s)
#     isos[i,2] =  vels[i].d_y/(u.km/u.s)
#     isos[i,3] =  vels[i].d_z/(u.km/u.s)
#     v = np.sqrt(isos[i,1]**2+isos[i,2]**2+isos[i,3]**2)
#     isos[i,4] =  v*(1+1774/(v**2))*(10**stars[i,8])/ESFquery(1000/stars[i,2], stars[i,8], stars[i,9])
#     # i+=1
#     #now with gravitational focussing
    
# # assert i == stars.shape[0]
# with open('/home/hopkinsl/GAIA/data/isos200.pickle', 'wb') as f:
#     pickle.dump(isos,f)

with open('/home/hopkinsl/GAIA/data/isos200.pickle', 'rb') as f:
    isos = pickle.load(f)

#Based on 2k2 solution in Bailer-Jones et al. (2018)
oumuamua = coord.SkyCoord(ra = 279.4752*u.deg,
                             dec = 33.8595*u.deg,
                             distance = 1*u.pc,
                             pm_ra_cosdec = 0*u.mas/u.yr,
                             pm_dec = 0*u.mas/u.yr,
                             radial_velocity = -26.420410*u.km/u.s
                             )
oumuamua_v = oumuamua.transform_to('galactic').velocity

# Based on Bailer-Jones et al. (2020)
borisov = coord.SkyCoord(ra = 32.79720*u.deg,
                             dec = 59.44014*u.deg,
                             distance = 1*u.pc,
                             pm_ra_cosdec = 0*u.mas/u.yr,
                             pm_dec = 0*u.mas/u.yr,
                             radial_velocity = -32.286894*u.km/u.s
                             )
borisov_v = borisov.transform_to('galactic').velocity

# Based on Schoenrich 2010
LSR = coord.SkyCoord(l = 0*u.deg,
                             b = 0*u.deg,
                             distance = 1*u.pc,
                             pm_l_cosb = 0*u.mas/u.yr,
                             pm_b = 0*u.mas/u.yr,
                             radial_velocity = 0*u.km/u.s,
                             frame='galacticlsr')
LSR_v = LSR.transform_to('galactic').velocity

#plots -neaten!
cbar = 'grey'
cmap2d = mpl.colormaps['binary']
cou = 'hotpink'
cbo = 'lime'
clsr = 'gold'
binArrays = [np.linspace(-100, 100, 100),np.linspace(-100, 50, 100)]
zbinArrays = [np.linspace(-100, 50, 100),np.linspace(-50, 50, 100)]



# fig, ax = plt.subplots()
# ax.hist2d(vels.d_x/(u.km/u.s), vels.d_y/(u.km/u.s), bins=binArrays, cmap=cmap)
# ax.set_xlim((-100, 100))
# ax.set_ylim((-100, 100))#unweighted stars


fig, ax = plt.subplots()
ax.hist(stars[:,8], bins=np.linspace(-1.5,0.7, 30), weights = starWeights/starWeights.sum(), facecolor=cbar)
ax.set_xlabel(r'$\mathrm{[M/H]}$')
ax.axvline(-0.4, c='turquoise')
ax.axvline( 0.4, c='turquoise')
# ax.set_ylabel('V')
ax.set_title('Metallicity Distribution Function')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/MH.png', dpi=300)

fig, ax = plt.subplots()
ax.hist(stars[:,9], bins=np.linspace(-1,15, 50), weights = starWeights/starWeights.sum(), facecolor=cbar)
ax.set_xlabel(r'$\mathrm{age}/\mathrm{Gyr}$')
# ax.set_ylabel('V')
ax.set_title('Stellar Age Distribution Function')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/starage.png', dpi=300)

fig, ax = plt.subplots()
ax.hist(isos[:,0], bins=np.linspace(0.0,0.6,20), weights = isos[:,4]/isos[:,4].sum(), facecolor=cbar)
ax.axvline(0.3, color=cbo, label='2I')
ax.set_xlabel(r'$f_\mathrm{H_2O}$')
# ax.set_ylabel('V')
ax.legend()
ax.set_title('ISO Composition Distribution')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/fH2O.png', dpi=300)

fig, ax = plt.subplots()
ax.hist2d(vels.d_x/(u.km/u.s), vels.d_y/(u.km/u.s), weights=starWeights, bins=binArrays, cmap=cmap2d)
ax.scatter(LSR_v.d_x, LSR_v.d_y, marker='^', color=clsr, label='LSR')
ax.set_xlabel(r'$U/\mathrm{kms}^{-1}$')
ax.set_ylabel(r'$V/\mathrm{kms}^{-1}$')
ax.legend()
ax.set_title('In-Plane Stellar Velocities')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/sUV.png', dpi=300)

fig, ax = plt.subplots()
ax.hist2d(vels.d_y/(u.km/u.s), vels.d_z/(u.km/u.s), weights=starWeights, bins=zbinArrays, cmap=cmap2d)
ax.scatter(LSR_v.d_y, LSR_v.d_z, marker='^', color=clsr, label='LSR')
ax.set_xlabel(r'$V/\mathrm{kms}^{-1}$')
ax.set_ylabel(r'$W/\mathrm{kms}^{-1}$')
ax.legend()
ax.set_title('Out-of-Plane Stellar Velocities')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/sVW.png', dpi=300)

fig, ax = plt.subplots()
ax.hist2d(isos[:,1], isos[:,2], bins=binArrays, weights = isos[:,4], cmap=cmap2d)
ax.scatter(oumuamua_v.d_x, oumuamua_v.d_y, marker='^', color=cou, label="1I")
ax.scatter(borisov_v.d_x, borisov_v.d_y, marker='^', color=cbo, label="2I")
ax.scatter(LSR_v.d_x, LSR_v.d_y, marker='^', color=clsr, label='LSR')
ax.set_xlabel(r'$U/\mathrm{kms}^{-1}$')
ax.set_ylabel(r'$V/\mathrm{kms}^{-1}$')
ax.legend()
ax.set_title('In-Plane ISO Velocities')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/iUV.png', dpi=300)

fig, ax = plt.subplots()
ax.hist2d(isos[:,2], isos[:,3], bins=zbinArrays, weights = isos[:,4], cmap=cmap2d)
ax.scatter(oumuamua_v.d_y, oumuamua_v.d_z, marker='^', color=cou, label='1I')
ax.scatter(borisov_v.d_y, borisov_v.d_z, marker='^', color=cbo, label='2I')
ax.scatter(LSR_v.d_y, LSR_v.d_z, marker='^', color=clsr, label='LSR')
ax.set_xlabel(r'$V/\mathrm{kms}^{-1}$')
ax.set_ylabel(r'$W/\mathrm{kms}^{-1}$')
ax.legend()
ax.set_title('Out-of-Plane ISO Velocities')
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/iVW.png', dpi=300)

#on-sky
def spherical(U,V,W):
    absv = np.sqrt(U**2 + V**2 + W**2)
    theta = np.arcsin(W/absv)
    phi = np.arctan2(V, U)
    return (absv, theta, phi)
speed, theta, phi = spherical(isos[:,1],isos[:,2],isos[:,3])

hporder=4

skybin_each_star = hp.ang2pix(hp.order2nside(hporder),180+phi*180/np.pi,-theta*180/np.pi, lonlat=True)
sky_weighted_counts = np.zeros(hp.order2npix(hporder))
for i in range(hp.order2npix(hporder)):
    sky_weighted_counts[i] = np.sum(isos[skybin_each_star==i,4])

# hp.mollview(sky_weighted_counts, cmap=cmap2d,cbar=False,title='ISO Radiants')
# hp.graticule( )

hp.newvisufunc.projview(
    sky_weighted_counts, graticule=True, graticule_labels=True,
    projection_type="mollweide", cmap=cmap2d, cbar=False,
    # custom_xtick_labels=['300°','240°','180°','120°','60°'], - Double check directions!
    xlabel=r"$l$",
    ylabel=r"$b$",
    title='ISO Radiant'
)
fig = plt.gcf()
ax = fig.get_axes()[0]
ouspeed, outheta, ouphi = spherical(oumuamua_v.d_x.value, oumuamua_v.d_y.value, oumuamua_v.d_z.value)
bospeed, botheta, bophi = spherical(borisov_v.d_x.value, borisov_v.d_y.value, borisov_v.d_z.value)
lsrspeed, lsrtheta, lsrphi = spherical(LSR_v.d_x.value, LSR_v.d_y.value, LSR_v.d_z.value)
# ax.scatter(-(np.pi+ouphi), -outheta, marker='^', s=100, color=cou, label='1I') #here - in phi is because ax object doesn't realise x-axis is backwards
# ax.scatter(-(np.pi+bophi), -botheta, marker='^', s=100, color=cbo, label='2I')
ax.scatter(-(np.pi+lsrphi), -lsrtheta, marker='^', s=100, color=clsr, label='Solar apex')
# ax.legend() causes crash
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/radiant_sans.png', dpi=300)



#chemodynamics - getting thirding
args = np.argsort(isos[:,0])
csum = np.cumsum(isos[args,4])
csum/=csum[-1]
masks = [args[csum<0.33], args[csum>0.66]]
names = ['Low','High']

for i in range(2): 
    m = masks[i]
    n = names[i]
    fig, ax = plt.subplots()
    ax.hist2d(isos[m,1], isos[m,2], bins=binArrays, weights = isos[m,4], cmap=cmap2d)
    ax.scatter(oumuamua_v.d_x, oumuamua_v.d_y, marker='^', color=cou, label='1I')
    ax.scatter(borisov_v.d_x, borisov_v.d_y, marker='^', color=cbo, label='2I')
    ax.scatter(LSR_v.d_x, LSR_v.d_y, marker='^', color=clsr, label='LSR')
    ax.set_xlabel(r'$U/\mathrm{kms}^{-1}$')
    ax.set_ylabel(r'$V/\mathrm{kms}^{-1}$')
    ax.set_title(f'In-Plane ISO Velocities - {n}')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'/home/hopkinsl/GAIA/plots/iUV_{n}.png', dpi=300)
    
    fig, ax = plt.subplots()
    ax.hist2d(isos[m,2], isos[m,3], bins=zbinArrays, weights = isos[m,4], cmap=cmap2d)
    ax.scatter(oumuamua_v.d_y, oumuamua_v.d_z, marker='^', color=cou, label='1I')
    ax.scatter(borisov_v.d_y, borisov_v.d_z, marker='^', color=cbo, label='2I')
    ax.scatter(LSR_v.d_y, LSR_v.d_z, marker='^', color=clsr, label='LSR')
    ax.set_xlabel(r'$V/\mathrm{kms}^{-1}$')
    ax.set_ylabel(r'$W/\mathrm{kms}^{-1}$')
    ax.set_title(f'Out-of-Plane ISO Velocities - {n}')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'/home/hopkinsl/GAIA/plots/iVW_{n}.png', dpi=300)
    
    skybin_each_star = hp.ang2pix(hp.order2nside(hporder),180+phi[m]*180/np.pi,-theta[m]*180/np.pi, lonlat=True)
    sky_weighted_counts = np.zeros(hp.order2npix(hporder))
    for i in range(hp.order2npix(hporder)):
        sky_weighted_counts[i] = np.sum(isos[m&(skybin_each_star==i),4])
    hp.newvisufunc.projview(
        sky_weighted_counts, graticule=True, graticule_labels=True,
        projection_type="mollweide", cmap=cmap2d, cbar=False,
        # custom_xtick_labels=['300°','240°','180°','120°','60°'], - Double check directions!
        xlabel=r"$l$",
        ylabel=r"$b$",
        title=f'ISO Radiant - {n} '+r'$f_\mathrm{H_2O}$'
    )
    fig = plt.gcf()
    ax = fig.get_axes()[0]
    ax.scatter(-(np.pi+ouphi), -outheta, marker='^', s=100, color=cou, label='1I')
    ax.scatter(-(np.pi+bophi), -botheta, marker='^', s=100, color=cbo, label='2I')
    ax.scatter(-(np.pi+lsrphi), -lsrtheta, marker='^', s=100, color=clsr, label='Solar apex')
    # ax.legend() causes crash
    fig.tight_layout()
    fig.savefig(f'/home/hopkinsl/GAIA/plots/radiant_{n}.png', dpi=300)


#Cool one - abs v dist
fig, ax = plt.subplots()
hist = ax.hist(np.sqrt(isos[:,1]**2+isos[:,2]**2+isos[:,3]**2), bins=np.linspace(0,150,51), weights = isos[:,4]/(3*isos[:,4].sum()), facecolor=cbar)
# ax.axvline(26.4, c=cou, label='1I')
# ax.axvline(32.2, c=cbo, label='2I')
# ax.axvline(0.8, c='royalblue', label='Arend–Roland')
ax.set_title('Predicted ISO speed distribution')
ax.set_xlabel(r'$v_\infty$ / $\mathrm{km}\mathrm{s}^{-1}$')
ax.legend()
fig.tight_layout()
fig.savefig('/home/hopkinsl/GAIA/plots/speed0.png', dpi=300)
# 3d v gives quite high range and may be unreliable due  reaching edges of disk - but this should reduce speeds?
#aybe weighting high speed thick disk stars too high?

# fig, ax = plt.subplots()
# ax.hist(np.sqrt(isos[:,1]**2+isos[:,2]**2), bins=np.linspace(0,150), weights = isos[:,4])
# # in-plane v is just as extended




















