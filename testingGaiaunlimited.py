import os

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt

from astroquery.gaia import Gaia

from gaiaunlimited import fetch_utils, utils, DR3SelectionFunctionTCG, subsample


from gaiaunlimited.utils import get_healpix_centers
import healpy as hp
#order these

assert Gaia.MAIN_GAIA_TABLE=='gaiadr3.gaia_source'

os.environ['SSL_CERT_FILE'] = '/Users/hopkinsm/GAIA/gaiaEnv/lib/python3.10/site-packages/certifi/cacert.pem'
# line needed for querying on mac

# mapHpx7 = DR3SelectionFunctionTCG()
# loads empirical, M10 based selection function from TCG et al. 2023
# options include hpx7 (default) which used order 7 healpix everywhere,
# 'multi' which use up to order 10 in dense regions, and 'patch' for 
# order 12 maps on the fly
# use is in .query(coords,G) method, which evaluates TCG SF at magnitude G and coord
# just make sure coords have at least one point per healpix 

# coords_of_centers = get_healpix_centers(7)
# gmag = np.ones_like(coords_of_centers) * 21.
# print(f'Computing the completeness for {len(coords_of_centers)} points.')
# completeness = mapHpx7.query(coords_of_centers,gmag)
# print(completeness)


# hp.mollview(completeness,
#             coord=['Celestial','Galactic'] # I think converts healpix index (corresponding to
#             # order in equatorial) to Galactic coords
#            )
# hp.mollview(completeness,
#             coord=['C','G'] # same as above. Due to confusion in labels, use above
#            )
# hp.mollview(completeness,
#             # shows default is Celestial (RA,dec), this is tilt of solar system
#            )
# hp.mollview(completeness,
#             coord=['C','E'] # mapping to Equatorial shows symmettry of scanning relative to Solar System
#            )
# hp.mollview(completeness,
#             coord=['G'] # shows putting one only changes label, doesn't do transform
#            )

def plotCompleteness(G, coords=get_healpix_centers(7), map=DR3SelectionFunctionTCG()):
    gArray = np.ones_like(coords)*G
    completeness = map.query(coords,gArray)
    hp.mollview(completeness, title=f"parent, G={G}",
                coord=['Celestial','Galactic']
               )
    
# plotCompleteness(21)
# GC fields remain at low completion even for magnitudes of 19 and less.
# Other than this sky complete at G=20.
# Between G=20 and G=21.5, completeness drops to essentially 0.
# RVS is limited to more like G<15, and GSP-Spec to G<14 

# RVS:
    
inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-1.3,5.1,0.4],
            'parallax': [10,30,10], # for actual shouldn't bin in parallax, for now this has effect of reducing time
         'mh_gspspec': [-0.9,0.9, 0.9] # for actual shouldn't bin here, just testing table join of ap worked, and if ap. or s. must be specified
        }

print('starting')
subsampleSF = subsample.SubsampleSelectionFunction(
                subsample_query = "radial_velocity is not null",
                file_name = 'dr3_lvl4', hplevel_and_binning = inDict)




