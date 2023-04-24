import os
import time
import numpy as np
from numpy import ma 
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
# from astropy.table import Table

import astropy.units as u
from astropy.coordinates import SkyCoord

<<<<<<< HEAD
os.environ['SSL_CERT_FILE'] = '/Users/hopkinsm/GAIA/gaiaEnv/lib/python3.10/site-packages/certifi/cacert.pem'
# line needed for querying
=======
# os.environ['SSL_CERT_FILE'] = '/Users/hopkinsm/GAIA/gaiaEnv/lib/python3.10/site-packages/certifi/cacert.pem'
# line needed for querying on mac
>>>>>>> b857449c9ffa48af2f3ec2b1f1c672a797ee4579

# coord = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
# width = u.Quantity(0.1, u.deg)
# height = u.Quantity(0.1, u.deg)
# r = Gaia.query_object_async(coordinate=coord, width=width, height=height)



tableNames = Gaia.load_tables(only_names=True)
for t in tableNames:
    if 'gcns' in t.get_qualified_name():
     # print(t.get_qualified_name())
     pass
 

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
