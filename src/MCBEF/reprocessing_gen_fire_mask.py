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

import pylib.MCBEF.MCBEF as MC
import numpy as np
import time
import xarray as xr

#---------------
# initialization 
#---------------
# initialize the namelist for algorithm configuration
nl = MC.NL.namelist_init('./namelist.input')

# initialized the time
# create the necessary time string needed in the detection
tt =  MC.ST.init_time(jdn, overpass_beg, overpass_end)

#-----------------------------
# preparing the satellite data
#-----------------------------
# read fire detections from the Level 2 dataset
filda_dict = MC.IO.read_fire_mask(nl, tt)


#-----------------------------
# write the data into netCDF
#-----------------------------
# Create an xarray Dataset
ds = xr.Dataset()
encodings = {}

dim_0 = {}
dim_1 = {}

# Define encoding options for compression and data type
encoding = {
	'zlib': True,          # Enable zlib compression
	'complevel': 9,        # Compression level (1-9)
	'dtype': 'int8'        # Ensuring data type is int8
}

attrs_legend = '\n0 not-processed (non-zero QF)\n1 bowtie\n2 glint\n3 water\n4 clouds\n5 clear land\n6 unclassified fire pixel\n7 low confidence fire pixel\n8 nominal confidence fire pixel\n9 high confidence fire pixel\n'

for key in filda_dict.keys():

	dim0, dim1 = filda_dict[key].shape
	if dim0 not in dim_0.keys():
		dim_0[dim0] = f'y{dim0}'
		
	if dim1 not in dim_1.keys():
		dim_1[dim1] = f'x{dim1}'
	
	# Adding DataArray to Dataset without defining coords
	ds[key] = xr.DataArray(filda_dict[key].astype(np.int8), dims=(dim_0[dim0], dim_1[dim1]))	
	ds[key].attrs['legend'] = attrs_legend
	
	
	encodings[key] = encoding
	

SAVE_PATH = '/Dedicated/jwang-data2/mzhou/project/OPNL_FILDA/BeamFire/REPROCESSING/FIRE_MASK/'
netcdf_path = SAVE_PATH + 'Fire_mask.' + sat + '.' + jdn + '.' + overpass_beg + '.' + overpass_end + '.nc'
# Write the dataset to a NetCDF file
ds.to_netcdf(path=netcdf_path, mode='w', encoding=encodings)












