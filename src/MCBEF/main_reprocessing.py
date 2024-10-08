import os
import sys
print(' - System: a python script has been called with paramters:', sys.argv)

# jdn          = 'A2019213'
# overpass_beg = '0842'
# overpass_end = '0849'
# sat 		 = 'VNP'

jdn          = sys.argv[1]
overpass_beg = sys.argv[2]
overpass_end = sys.argv[3]
sat 		 = sys.argv[4]

# initiate the compiling path...
os.environ['THEANO_FLAGS'] = f'base_compiledir=./THEANO/theano_{jdn}'
os.environ['PYTENSOR_FLAGS'] = f'base_compiledir=./PYTENSOR/pytensor_{jdn}'

import pylib.MCBEF.MCBEF as MC
import numpy as np
import time


# initialize the namelist for algorithm configuration
nl = MC.NL.namelist_init('./namelist.input')

# initialized the time
# create the necessary time string needed in the detection
tt =  MC.ST.init_time(jdn, overpass_beg, overpass_end)

#-----------------------------
# preparing the satellite data
#-----------------------------
filda_dict = MC.IO.read_data_reprocess(nl, tt)

#---------------------
# preparing PYMC model
#---------------------
# initialize the sensor
v_sensor = MC.SR.viirs()

# define the background sensor (bgs)
bgs = MC.SR.GetSensor(v_sensor, nl.sel_bg_bands)

# define the fire sensitive sensor (fss)
fss = MC.SR.GetSensor(v_sensor, nl.sel_fire_bands)

# initialize the background temperature estimator
bg_estimator, ts_bg = MC.init_background_estimation(nl, bgs)

# initialize the background temperature estimator
biphase_estimator, ts_bi_fire = MC.init_biphasic_estimation(nl, fss)

# initialize the background temperature estimator
uniphase_estimator, ts_uni_fire = MC.init_uniphasic_estimation(nl, fss)

#---------
# sampling
#---------
# define two dataset object for data processing
data_fire, data_bg = MC.get_sample_set(filda_dict, nl)

print(' - MAIN: Number of fire', len(data_fire))
start_time_main = time.time()
result = MC.MCBEF_MP(nl, data_bg, data_fire, 
                     bg_estimator, biphase_estimator, uniphase_estimator,
                     ts_bg, ts_bi_fire, ts_uni_fire,
                     bgs, fss)
end_time_main = time.time()  # End time measurement
execution_time_main = end_time_main - start_time_main  # Calculate execution time
print(' - MAIN: Number of estimation', len(result))
print(f" - MAIN: Execution time {execution_time_main:4.2f} seconds\n")

#----------------
# post-processing 
#----------------
sample_dicts, state_dict = MC.post_processing(filda_dict, result, nl)

#-----------------
# write the output
#-----------------
sample_path, state_path = MC.IO.ini_output_dir(nl, tt)

for op in sample_dicts.keys():
    savename= sample_path + sat + '47MCBEF.Sample.' + jdn + '.' + \
              op + '.' + nl.version + '.nc'
    MC.IO.dict2nc(sample_dicts[op], nl, sat, 'Sample', savename)

savename=state_path + sat + '47MCBEF.State.' + jdn + '.' + \
         overpass_beg+'_'+overpass_end + '.' + nl.version + '.nc'
MC.IO.dict2nc(state_dict, nl, sat, 'State', savename)


