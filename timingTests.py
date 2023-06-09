import time
import os
import numpy as np
from numpy import ma 
from astroquery.gaia import Gaia


"""To test
source vs soure_lite
join vs no join
finding one bin by using WHERE v finding multiple bins using floor and indexes

for all for now redcing number of bins to reduce time

also test effectd of single quotes in not null

Results
quotes needed! gues its to do with evaluation in i_t_e func

in a trial of one reea, source_lite was three times slower than source. may be due to dodgy train wifi
then got 104 error connection reset by peer? very quikyl approx 5 sec
source is on 45, 53
source_lite i on 136, 37

got one time
IncompleteRead: IncompleteRead(556512 bytes read, 1809323 more expected)
something goeing wrog, either mem leak or bad train wii, I'll try again later

Likely possibilities
2<1, 4<3, 6<5 - source_lite is just faster than source
4<3, 6<5 - source_lite is faster than source when joining
6<5 - source_lite is faster than source when joining and using
Jump from 1,2 to 3,4> jump from 3,4 to 5,6 means join is slow
Jump from 1,2 to 3,4< jump from 3,4 to 5,6 means EITHER calculations are slow
OR join is slow but clever, so doesn't occur when not needed
A sign of this would be 1==3 and 2==4, and  in the times of 7-10, which dont do any calculations
but use/dont use ligte and d/dont use mh

Results:
    first querys are always slower - more representative? what makes later onesfaster?
unclear but I don't think it's cached locally
restanting computer doesn't reset it, later queries are also fast and 
changing query slightly doesn't reset it so there must be some element saved o make it fast

Alternatively, give up trying to understand and just try the query I want one bin at a time'

"tuns out, can't get precaculated colour in source_lite"
"""

temp0 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source_lite WHERE phot_g_mean_mag<12 and parallax>5"""



query1 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source WHERE phot_g_mean_mag<12 and parallax>5"""
#defalt

query2 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source_lite WHERE phot_g_mean_mag<12 and parallax>5.1"""
#source_lite

query3 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
 WHERE phot_g_mean_mag<12 and parallax>5"""
#join but dont use and in selection

query4 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL',1.0,0.0))
FROM gaiadr3.gaia_source_lite JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE phot_g_mean_mag<12 and parallax>5"""
#join with lite but dont use and in selection      

query5 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
 WHERE phot_g_mean_mag<12 and parallax>5"""
#join, use and in selection

query6 = """SELECT source_id, to_integer(IF_THEN_ELSE('radial_velocity IS NOT NULL AND mh_gspspec IS NOT NULL',1.0,0.0))
 FROM gaiadr3.gaia_source_lite JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE phot_g_mean_mag<12 and parallax>5"""
#join with lite, use and in selection     

query7 = """SELECT source_id, radial_velocity
 FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
 WHERE phot_g_mean_mag<12 and parallax>5"""
#join but dont use and in selection and doesnt do any calculations

query8 = """SELECT source_id, radial_velocity
 FROM gaiadr3.gaia_source_lite JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE phot_g_mean_mag<12 and parallax>5"""
#join with lite but dont use and doesnt do any calculations

query9 = """SELECT source_id, radial_velocity, mh_gspspec
 FROM gaiadr3.gaia_source JOIN gaiadr3.astrophysical_parameters USING (source_id)
 WHERE phot_g_mean_mag<12 and parallax>5"""
#join but does use and in selection and doesnt do any calculations

query10 = """SELECT source_id, radial_velocity, mh_gspspec
 FROM gaiadr3.gaia_source_lite  JOIN gaiadr3.astrophysical_parameters USING (source_id)
WHERE phot_g_mean_mag<12 and parallax>5"""
#join with lite but does use and doesnt do any calculations


query = query9
print(query)

t0 = time.time()
ts = [t0]
for i in range(1):
    job = Gaia.launch_job_async(query, verbose=True)
    ts.append(time.time())
    
    table = job.get_results()
    print('N = ', len(table))
    print(f'time = {ts[-1]-ts[-2]} seconds')