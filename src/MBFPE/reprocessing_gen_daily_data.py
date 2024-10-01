import os
import sys
print(' - System: a python script has been called with paramters:', sys.argv)

# jdn          = 'A2019213'
# overpass_beg = '0842'
# overpass_end = '0849'
# sat 		 = 'VNP'
# python reprocessing_gen_daily_data.py A2019213 0842 0849 VNP

jdn          = sys.argv[1]
overpass_beg = sys.argv[2]
overpass_end = sys.argv[3]
sat 		 = sys.argv[4]

import pylib.MBFPE.MBFPE as MB
import numpy as np
import time
import xarray as xr

#---------------
# initialization 
#---------------
# initialize the namelist for algorithm configuration
nl = MB.NL.namelist_init('./namelist.input')

# initialized the time
# create the necessary time string needed in the detection
tt =  MB.ST.init_time(jdn, overpass_beg, overpass_end)

#-----------------------------
# preparing the satellite data
#-----------------------------
# read fire detections from the Level 2 dataset
filda_dict = MB.IO.read_data(nl, tt)

# add the surface emissivity
filda_dict = MB.IO.get_surface_emit(filda_dict, nl, tt)

# add the static flag
filda_dict = MB.IO.get_static_thermal_anomaly(filda_dict, nl, tt)

# add the land surface type
filda_dict = MB.IO.get_surface_type(filda_dict, nl, tt)

print(filda_dict.shape)
filda_dict.index = np.arange(filda_dict.shape[0])


SAVE_PATH = '/Dedicated/jwang-data2/mzhou/project/OPNL_FILDA/BeamFire/REPROCESSING/DATA/'
# Convert DataFrame to xarray Dataset
ds = xr.Dataset.from_dataframe(filda_dict)
# Write the dataset to a NetCDF file
netcdf_path = SAVE_PATH + 'Intermediate.' + sat + '.' + jdn + '.' + overpass_beg + '.' + overpass_end + '.nc'
ds.to_netcdf(path=netcdf_path, mode='w')
print("NetCDF file has been created:", netcdf_path)
