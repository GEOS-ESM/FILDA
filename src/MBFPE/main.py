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
# os.environ['THEANO_FLAGS'] = 'base_compiledir=./THEANO/theano_{}'.format('A2019220')
os.environ['THEANO_FLAGS'] = f'base_compiledir=./THEANO/theano_{jdn}'
os.environ['PYTENSOR_FLAGS'] = f'base_compiledir=./PYTENSOR/pytensor_{jdn}'

import pylib.MBFPE.MBFPE as MB
import numpy as np
import time

#---------------
# initialization 
#---------------
# initialize the namelist for algorithm configuration
nl = MB.NL.namelist_init('./namelist.input')

# initialized the time
# create the necessary time string needed in the detection
tt =  MB.ST.init_time(jdn, overpass_beg, overpass_end)

#---------------------
# preparing PYMC model
#---------------------
# initialize the sensor
v_sensor = MB.SR.viirs()

# define the background sensor (bgs)
bgs = MB.SR.GetSensor(v_sensor, nl.sel_bg_bands)

# define the fire sensitive sensor (fss)
fss = MB.SR.GetSensor(v_sensor, nl.sel_fire_bands)

# initialize the background temperature estimator
bg_estimator, ts_bg = MB.init_background_estimation(nl, bgs)

# initialize the background temperature estimator
biphase_estimator, ts_bi_fire = MB.init_biphasic_estimation(nl, fss)

# initialize the background temperature estimator
uniphase_estimator, ts_uni_fire = MB.init_uniphasic_estimation(nl, fss)

#-----------------------------
# preparing the satellite data
#-----------------------------
# read fire detections from the Level 2 dataset
filda_dict = MB.IO.read_data(nl, tt)

# add the surface emissivity
filda_dict = MB.IO.get_surface_emit(filda_dict, nl, tt)

# add the static flag
filda_dict = MB.IO.get_static_thermal_anomaly(filda_dict, nl, tt)

#---------
# sampling
#---------
# define two dataset object for data processing
data_fire, data_bg = MB.get_sample_set(filda_dict, nl)

print(' - MAIN: Number of fire', len(data_fire))
start_time_main = time.time()
result = MB.MBFFP_MP(nl, data_bg, data_fire, 
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
sample_dicts, state_dict = MB.post_processing(filda_dict, result, nl)

#-----------------
# write the output
#-----------------
sample_path, state_path = MB.IO.ini_output_dir(nl, tt)

for op in sample_dicts.keys():
    savename= sample_path + sat + '.Sample.' + jdn + '.' + \
              op + '.' + nl.version + '.nc'
    MB.IO.dict2nc(sample_dicts[op], nl, sat, 'Sample', savename)

savename=state_path + sat + '.State.' + jdn + '.' + \
         overpass_beg+'_'+overpass_end + '.' + nl.version + '.nc'
MB.IO.dict2nc(state_dict, nl, sat, 'State', savename)







