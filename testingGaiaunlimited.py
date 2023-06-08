import os
import time

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
    
# inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-1.3,5.1,0.4],#[-2.5,5.1,0.4],
#             'parallax': [10,30,10], # for actual shouldn't bin in parallax, for now this has effect of reducing time
#          'mh_gspspec': [-0.1,0.1, 0.1] # for actual shouldn't bin here, just testing table join of ap worked, and if ap. or s. must be specified
#         }
# fname = 'test'

# inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-1.3,5.1,0.4],#[-2.5,5.1,0.4],
#             'parallax': [10,30,10], # for actual shouldn't bin in parallax, for now this has effect of reducing time
#          'mh_gspspec': [-0.1,0.1, 0.1] # for actual shouldn't bin here, just testing table join of ap worked, and if ap. or s. must be specified
#         }
# fname = 'hp5, '

# inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-1.3,5.1,0.4],#[-2.5,5.1,0.4],
#             'parallax': [10,30,10], # for actual shouldn't bin in parallax, for now this has effect of reducing time
#          'mh_gspspec': [-0.1,0.1, 0.1] # for actual shouldn't bin here, just testing table join of ap worked, and if ap. or s. must be specified
#         }
# fname = 'test'

# print('starting')
# subsampleSF = subsample.SubsampleSelectionFunction(
#                 subsample_query = "radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL",
#                 file_name = 'lvl4', hplevel_and_binning = inDict, use_astrophysical_parameters=True)

#I'll run this with SF in paper to check, and one with mh is not null,
#and one using builti rv selection(also hp5)
# and finally one which should cause error
#name well, and wll be saved locally
# would be nice to have fig 5 type plot - unclear how they made it, but 
# probalby can make by n-weighting healpix
#bug shouldn't cause problem if I follow paper - therefore lvl4

startTime=time.time()

try:
    
    # print('run1')
    # inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-2.5,5.1,0.4]}
    # fname = 'hp4_docsRange_rv'
    # ssquery = "radial_velocity IS NOT NULL"
    # subsampleSF = subsample.SubsampleSelectionFunction(
    #                 subsample_query = ssquery,
    #                 file_name = fname, hplevel_and_binning = inDict,
    #                 use_astrophysical_parameters=False)
    
    print('run2')
    inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-2.5,5.1,0.4]}
    fname = 'hp4_docsRange_rv_mh'
    ssquery = "radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL"
    subsampleSF = subsample.SubsampleSelectionFunction(
                    subsample_query = ssquery,
                    file_name = fname, hplevel_and_binning = inDict,
                    use_astrophysical_parameters=True)
    
    # print('run3')
    # rvssf = subsample.DR3RVSSelectionFunction()
    
    # print('run4')
    # inDict = {'healpix': 4,'phot_g_mean_mag': [3,20,0.2],'g_rp': [-2.5,5.1,0.4]}
    # fname = 'hp4_docsRange_rv_mh_noAP'
    # ssquery = "radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL"
    # subsampleSF = subsample.SubsampleSelectionFunction(
    #                 subsample_query = ssquery,
    #                 file_name = fname, hplevel_and_binning = inDict,
    #                 use_astrophysical_parameters=False)
    
except ConnectionResetError:
    print(time.time()-startTime)

print(time.time()-startTime)


#explaination of how query works
"""once filled in, it's 
SELECT healpix_, phot_g_mean_mag_, g_rp_, COUNT(*) AS n, SUM(selection) AS k
FROM (
      SELECT 
            to_integer(GAIA_HEALPIX_INDEX(4,source_id)) AS healpix_,
            to_integer(floor((phot_g_mean_mag - 3)/0.2)) AS phot_g_mean_mag_,
            to_integer(floor((g_rp - -1.3)/0.4)) AS g_rp_,
            to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL 
                                    AND mh_gspspec IS NOT NULL',1.0,0.0)) AS selection
      FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
      WHERE phot_g_mean_mag > 3    AND phot_g_mean_mag < 20 
            AND        g_rp > -1.3 AND            g_rp < 5.1
     ) AS subquery
GROUP BY healpix_, phot_g_mean_mag_, g_rp_

Firstly, FROM (...) AS subquery defines a subquery which builds a
'temporary table' that the outer clauses draw from.
In subquery, goes through gaiadr3.source (and ap), limited to stars in the 
binning range, recording the bin index of each star (floor etc)
as healpix_, phot_g_mean_mag_, and g_rp_, and recording if rv and mh are not 
null in a column called 'selection'.
 Result is that the 'temporary table' has row for each star in binning
range, with the bin index in healpix_, phot_g_mean_mag_, g_rp_, and if in subsample in selection
Outside subquery, GROUP BY means that for each combo value of  healpix_, phot_g_mean_mag_, g_rp_
(i.e bin) counts number of rows as n and number in subsample (selection=1) as k

This is why if no stars in bin it isn' carried forward. There aren't the rows to group 
"""

"""
looks like everything working. one tha meat to fail did, and precompted rvssf looks like my version 
slihgt differences attributed to mine being hp4, theirs hp5
mh version timed out - could it be becasue of join to ap?
afer how long did it time out? if current run still going at 11.15 its an hour
can it b sped up?
or broken up: could see where problems are

Running with timer: if exactly an hour (11.35 start) then likely a tieout problem
https://www.cosmos.esa.int/web/gaia-users/archive/faq
says anonymous tiout is 90min, registered and logged in is 120min

Depending on what is slow bit, could use source_lite and/or count(healpix_) instead of count(*)
assuming no NULL in healPix_

could possibly just do part of work in query, rest locally. Or all, on cluster

slow bit is probably subquery:
SELECT 
      to_integer(GAIA_HEALPIX_INDEX(4,source_id)) AS healpix_,
      to_integer(floor((phot_g_mean_mag - 3)/0.2)) AS phot_g_mean_mag_,
      to_integer(floor((g_rp - -1.3)/0.4)) AS g_rp_,
      to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL 
                              AND mh_gspspec IS NOT NULL',1.0,0.0)) AS selection
FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE phot_g_mean_mag > 3    AND phot_g_mean_mag < 20 
      AND        g_rp > -1.3 AND            g_rp < 5.1

by advice on gaia website calculaions in query are slow, could get all data and calc bins on cluster

OR do each bin as a seprate quesry! could deal with n=0 bug while I'm at it
at 100 minutes, hasn't timed out, but last timefailed at 60ish. Unsure what's going wrong, but splitting it up will help
"""
