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
    
"""
returning to original subsample query strategy as that is faster"""
    

"""Run 1 rv gave error
***Query 1520/1530 done, time elapsed = 43438.91185951233 sec***
Traceback (most recent call last):

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/spyder_kernels/py3compat.py:356 in compat_exec
    exec(code, globals, locals)

  File ~/GAIA/iso-gaia/makeSubsampleSelFunc.py:42
    subsampleSF = subsample.SubsampleSelectionFunction(

  File ~/GAIA/Packages/gaiaunlimited/src/gaiaunlimited/selectionfunctions/subsample.py:306 in __init__
    df = _download_binned_subset(self)

  File ~/GAIA/Packages/gaiaunlimited/src/gaiaunlimited/selectionfunctions/subsample.py:291 in _download_binned_subset
    df[columns].to_csv(f, index = False)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/pandas/core/frame.py:3767 in __getitem__
    indexer = self.columns._get_indexer_strict(key, "columns")[1]

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/pandas/core/indexes/base.py:5876 in _get_indexer_strict
    self._raise_if_missing(keyarr, indexer, axis_name)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/pandas/core/indexes/base.py:5938 in _raise_if_missing
    raise KeyError(f"{not_found} not in index")

KeyError: "['phot_g_mean_mag_', 'g_rp_'] not in index"




Run 2 gave error mh

***Query 1350/1530 done, time elapsed = 27866.89084625244 sec***
Traceback (most recent call last):

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/spyder_kernels/py3compat.py:356 in compat_exec
    exec(code, globals, locals)

  File ~/GAIA/iso-gaia/makeSubsampleSelFunc.py:43
    subsampleSF = subsample.SubsampleSelectionFunction(

  File ~/GAIA/Packages/gaiaunlimited/src/gaiaunlimited/selectionfunctions/subsample.py:306 in __init__
    df = _download_binned_subset(self)

  File ~/GAIA/Packages/gaiaunlimited/src/gaiaunlimited/selectionfunctions/subsample.py:283 in _download_binned_subset
    job = Gaia.launch_job_async(query_to_gaia, name=self.file_name)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/gaia/core.py:903 in launch_job_async
    return TapPlus.launch_job_async(self, query=query,

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/core.py:459 in launch_job_async
    job.get_results()

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/model/job.py:246 in get_results
    self.__load_async_job_results()

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/model/job.py:342 in __load_async_job_results
    wjResponse, phase = self.wait_for_job_end()

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/model/job.py:327 in wait_for_job_end
    responseData = self.get_phase(update=True)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/model/job.py:185 in get_phase
    response = self.connHandler.execute_tapget(phase_request)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/conn/tapconn.py:196 in execute_tapget
    return self.__execute_get(context, verbose)

  File ~/GAIA/gaiaEnv/lib/python3.8/site-packages/astroquery/utils/tap/conn/tapconn.py:243 in __execute_get
    conn.request("GET", context, None, self.__getHeaders)

  File /usr/lib/python3.8/http/client.py:1256 in request
    self._send_request(method, url, body, headers, encode_chunked)

  File /usr/lib/python3.8/http/client.py:1302 in _send_request
    self.endheaders(body, encode_chunked=encode_chunked)

  File /usr/lib/python3.8/http/client.py:1251 in endheaders
    self._send_output(message_body, encode_chunked=encode_chunked)

  File /usr/lib/python3.8/http/client.py:1011 in _send_output
    self.send(msg)

  File /usr/lib/python3.8/http/client.py:951 in send
    self.connect()

  File /usr/lib/python3.8/http/client.py:1425 in connect
    self.sock = self._context.wrap_socket(self.sock,

  File /usr/lib/python3.8/ssl.py:500 in wrap_socket
    return self.sslsocket_class._create(

  File /usr/lib/python3.8/ssl.py:1040 in _create
    self.do_handshake()

  File /usr/lib/python3.8/ssl.py:1309 in do_handshake
    self._sslobj.do_handshake()

TimeoutError: [Errno 110] Connection timed out
"""




