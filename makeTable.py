import time
import numpy as np
from numpy import ma 
from astroquery.gaia import Gaia
import pickle


query = ('SELECT ap.mh_gspspec AS mh, s.parallax AS plx, s.ra AS ra, s.dec AS dec, s.pmra AS pmra, s.pmdec AS pmdec,  s.radial_velocity AS rv '
          +'FROM gaiadr3.gaia_source_lite AS s '
          +'JOIN gaiadr3.astrophysical_parameters AS ap USING (source_id) '
          +'WHERE s.parallax>10 AND ap.mh_gspspec IS NOT NULL AND s.radial_velocity IS NOT NULL'
          ) # not sure how some stars have spectral mh but not rv, bu it happens quite frequently
print(query)

t1 = time.time()
job = Gaia.launch_job_async(query)
t2 = time.time()
print(f'time = {t2-t1}')
    
table = job.get_results()
print(table)
print(len(table))
print(table[1]['mh'] is ma.masked)

with open('data.pickle','wb') as f:
    pickle.dump(table, f)

# table.write('data.dat', format='ascii', overwrite=True)