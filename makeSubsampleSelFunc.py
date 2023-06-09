import os
import time
from math import isclose

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt

from astroquery.gaia import Gaia
import logging
logger = logging.getLogger("astroquery")
logger.setLevel(logging.WARNING)

from gaiaunlimited import fetch_utils, utils, DR3SelectionFunctionTCG, subsample

from gaiaunlimited.utils import get_healpix_centers


assert Gaia.MAIN_GAIA_TABLE=='gaiadr3.gaia_source'
os.environ['SSL_CERT_FILE'] = '/Users/hopkinsm/GAIA/gaiaEnv/lib/python3.10/site-packages/certifi/cacert.pem'
# line needed for querying on mac

# def arr(start, stop, step):
#     arr = np.arange(round((stop-start)/step)+1)*step+start
#     assert isclose(arr[-1], stop)
#     return arr

hplvl = 4
# gArr = arr(3,20,0.2)
# cArr = arr(-2.5,5.1,0.4)
# shape = ((len(gArr)-1), (len(cArr)-1))

fname = 'hp4_docsRange_rv_mh'
ssquery = "radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL"

# fname = 'hp4_docsRange_rv_split'
# ssquery = "radial_velocity IS NOT NULL"

inDict = {'healpix': hplvl,'phot_g_mean_mag': (3,20,0.2), 'g_rp': (-2.5,5.1,0.4)}
startTime=time.time()
try:
    subsampleSF = subsample.SubsampleSelectionFunction(
                    subsample_query = ssquery,
                    file_name = fname, hplevel_and_binning = inDict,
                    use_astrophysical_parameters=True)
    
except ConnectionResetError:
    print(time.time()-startTime)

print(time.time()-startTime)
    
    







