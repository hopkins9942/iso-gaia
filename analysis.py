import pickle
import numpy as np
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
# for now use 2 for plx, 4-6 for ra,dec,pmra, pmdec, 7 for ra, 8-9 for mh,age
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
#     isos[i,4] =  (10**stars[i,8])/ESFquery(1000/stars[i,2], stars[i,8], stars[i,9])
#     # i+=1
    
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

#plots -neaten!

fig, ax = plt.subplots()
ax.hist2d(vels.d_x/(u.km/u.s), vels.d_y/(u.km/u.s), bins=100)
ax.set_xlim((-100, 100))
ax.set_ylim((-100, 100))#unweighted stars


fig, ax = plt.subplots()
ax.hist2d(vels.d_x/(u.km/u.s), vels.d_y/(u.km/u.s), weights=starWeights, bins=100)
ax.set_xlim((-100, 100))
ax.set_ylim((-100, 100))#unweighted stars

fig, ax = plt.subplots()
ax.hist2d(isos[:,1], isos[:,2], bins=100)
ax.set_xlim((-100, 100))
ax.set_ylim((-100, 100))

fig, ax = plt.subplots()
ax.hist2d(isos[:,1], isos[:,2], bins=100, weights = isos[:,4])
ax.set_xlim((-100, 100))
ax.set_ylim((-100, 100)) # isos otably more concentratednear zero due to high FeH stars

fig, ax = plt.subplots()
ax.hist(vels.d_x/(u.km/u.s), bins=np.linspace(-100,100), weights = starWeights)

fig, ax = plt.subplots()
ax.hist(isos[:,1], bins=np.linspace(-100,100), weights = isos[:,4])

fig, ax = plt.subplots()
ax.hist(vels.d_z/(u.km/u.s), bins=np.linspace(-50,50), weights = starWeights)

fig, ax = plt.subplots()
ax.hist(isos[:,3], bins=np.linspace(-50,50), weights = isos[:,4])
#differences quite slight but that's fine

#Cool one - abs v dist

fig, ax = plt.subplots()
ax.hist(np.sqrt(isos[:,1]**2+isos[:,2]**2+isos[:,3]**2), bins=np.linspace(0,150), weights = isos[:,4])
# 3d v gives quite high range and may be unreliable due  reaching edges of disk - but this should reduce speeds?
#aybe weighting high speed thick disk stars too high?

fig, ax = plt.subplots()
ax.hist(np.sqrt(isos[:,1]**2+isos[:,2]**2), bins=np.linspace(0,150), weights = isos[:,4])
# in-plane v is just as extended




















