import os
import time
import numpy as np
from numpy import ma 
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
# from astropy.table import Table

import astropy.units as u
from astropy.coordinates import SkyCoord


os.environ['SSL_CERT_FILE'] = '/Users/hopkinsm/GAIA/gaiaEnv/lib/python3.10/site-packages/certifi/cacert.pem'
# line needed for querying on mac

# coord = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
# width = u.Quantity(0.1, u.deg)
# height = u.Quantity(0.1, u.deg)
# r = Gaia.query_object_async(coordinate=coord, width=width, height=height)



# tableNames = Gaia.load_tables(only_names=True)
# for t in tableNames:
#     if 'gcns' in t.get_qualified_name():
#         print(t.get_qualified_name())
#         # pass Confirmed gcns was in edr3 - weirdly couldn't find it last time
 

# job1 = Gaia.launch_job_async('SELECT source_id FROM gaiadr3.gaia_source_lite WHERE parallax>50')



# job = Gaia.launch_job_async('SELECT TOP 30 mh_gspspec FROM gaiadr3.astrophysical_parameters')
# print(job.get_results().colnames)
# table = job.get_results()
# print(table)


query = ('SELECT ap.mh_gspspec AS mh, s.parallax AS plx '
          +'FROM gaiadr3.gaia_source_lite AS s '
          +'JOIN gaiadr3.astrophysical_parameters AS ap USING (source_id) '
          +'WHERE s.parallax>10 AND ap.mh_gspspec IS NOT NULL'
          )
print(query)

t1 = time.time()
job = Gaia.launch_job_async(query)
t2 = time.time()
print(f'time = {t2-t1}')
    
table = job.get_results()
print(table)
print(len(table))
print(table[1]['mh'] is ma.masked)

print(f"number with mh vals: {table['mh'].count()}")

# output:

# time = 155.26039505004883
#   mh         plx        
#  dex         mas        
# ----- ------------------
#  0.32 10.000031865768985
# -0.32 10.000112983745307
# -0.11  10.00013429110453
#  0.02 10.000299262973787
#  0.18 10.000420315932535
# -0.19 10.000485084597802
#  0.09 10.000571114705792
# -0.86 10.000579426263918
#  0.18  10.00060517015851
# -1.32 10.000630749357805
#                ...
# -0.88  283.8401180184389
#  -0.5 285.99494829578117
# -0.52  286.0053518616485
# -1.02  296.3053079139394
#  -0.6  304.1353692001036
# -3.03 316.48118678226916
# -1.51  336.0266016683708
# -1.29 392.75294543876464
#  -3.0   546.975939730948
#  0.04  768.0665391873573
# Length = 81070 rows
# 81070
# False
# number with mh vals: 81070


