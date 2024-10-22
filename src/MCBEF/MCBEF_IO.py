'''

Python library of Multichannel Biphasic Fire Parameter (MCBEF)

IO module, v1.0


This module provide the basic reading and writing functions for MCBEF

'''

import numpy as np
import copy
import glob
from netCDF4 import Dataset
import sys, os
import time
import pandas as pd
from .MCBEF_utils import *
from .MCBEF_SPACE_TIME import *
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

def get_FILDA_size(filename, verbose=0):
	'''
	
	default setup to mute the message...
	
	'''	
	
	message = f'Reading {filename}...'
	printf(message, verbose, prefix = 'MCBEF_IO')
	
	ncid = Dataset(filename, 'r')
	num_fire = ncid.dimensions['nFire'].size
	
	ncid.close()
	
	return num_fire


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

def read_FILDA(filename, ignore_2d = True, verbose=0):
	'''
	'''
	import numpy as np
	from netCDF4 import Dataset	
	
	message = f'Reading {filename}...'
	printf(message, verbose, prefix = 'MCBEF_IO')
	
	output ={}
	
	ncid = Dataset(filename)
	params_lst = list(ncid.variables)

	if ignore_2d:
		for param in params_lst:
			n_dim = ncid[param][:].ndim
			if n_dim == 1:
				output[param] = ncid[param][:]
	else:
		for param in params_lst:
			output[param] = ncid[param][:]
			
	
	ncid.close()
	
	return output

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def gen_overpass(tt):

	minutes_beg = int(float(tt.OP1[0:2]) * 60 + float(tt.OP1[2:]))
	minutes_end = int(float(tt.OP2[0:2]) * 60 + float(tt.OP2[2:]))
	
	minutes = np.arange(minutes_beg, minutes_end, 6)
	
	ops = [ str(m//60).zfill(2)+str(m%60).zfill(2)  for m in minutes]
	
	return ops
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_data(nl, tt):
	
	ops = gen_overpass(tt)

	filda_dict = []
	
	for op in ops:

		file_temp = glob.glob(nl.input_path + 'VNP47MOD' + '/'+ tt.Y + \
		                      '/' + tt.DOY + '/' + '*' + op + '*.nc')
		                      
		if len(file_temp) > 0:
			message = f'Reading {tt.JDN} overpass {op}'
			printf(message, 1, prefix = 'MCBEF_IO')
			temp_dict = pd.DataFrame.from_dict(read_FILDA(file_temp[0]),'columns')
			temp_dict['overpass'] = op
			filda_dict.append(temp_dict)
			
		else:
			message = f'Cannot find {tt.JDN} overpass {op}'
			printf(message, 1, prefix = 'MCBEF_IO')
			
	if len(filda_dict) > 0:
		filda_dict = pd.concat(filda_dict)
		filda_dict['FP_DNB_Rad'] = filda_dict['DNB_observations']*1e-5
		filda_dict['FP_DNB_Rad_Mean'] = 0.0
		
	else:
		filda_dict = None	

	return filda_dict

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_data_reprocess(nl, tt):

    file_temp = glob.glob(nl.input_path + 'REPROCESSING/DATA/'+ '*' + \
                          tt.JDN + '*.nc')

    if len(file_temp) > 0:
        message = f'Reading {tt.JDN} reprocessing dataset'
        printf(message, 1, prefix = 'MCBEF_IO')
        
        filda_dict = pd.DataFrame.from_dict(read_FILDA(file_temp[0]),'columns')
        filda_dict.drop('index', axis=1, inplace=True)
    else:
        message = f'Cannot find {tt.JDN} reprocessing dataset'
        printf(message, 1, prefix = 'MCBEF_IO')
        filda_dict = None
        
    return filda_dict

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_data_viirs_dummy(nl):

	input_path = os.path.dirname(os.path.abspath(__file__))
	# Get the absolute path to the directory where the current 
	file_temp =  input_path + '/sensor/sensor_viirs_dummy_obs/VIIRS_obs_dummy.csv'
	
	filda_dict = pd.read_csv(file_temp)
	
	return filda_dict
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_fire_mask(nl, tt):
	
	ops = gen_overpass(tt)

	filda_dict = {}
	
	for op in ops:

		file_temp = glob.glob(nl.input_path + 'VNP47MOD' + '/'+ tt.Y + \
		                      '/' + tt.DOY + '/' + '*' + op + '*.nc')
		                      
		if len(file_temp) > 0:
			message = f'Reading {tt.JDN} overpass {op}'
			printf(message, 1, prefix = 'MCBEF_IO')
			temp_dict = pd.DataFrame.from_dict(read_FILDA(file_temp[0]),'columns')
			
			ncid = Dataset(file_temp[0], 'r')
			ncid.set_auto_mask(False)
			fire_mask = ncid['Fire_mask'][:]
			ncid.close()
			
			filda_dict[op] = fire_mask
			

		else:
			message = f'Cannot find {tt.JDN} overpass {op}'
			printf(message, 1, prefix = 'MCBEF_IO')

	return filda_dict
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_surface_emit(filda_dict, nl, tt, **kwargs):
	'''
	''' 
	numCeil	  = kwargs.get('numCeil', 1200)
	
	num_fire = filda_dict.shape[0]
	
	emit_params = ['Emis_M14', 'Emis_M15', 'Emis_M16', 'Emis_I05']
	
	# in any case, we need the climatology
	emis_dict = get_surface_emit_sinu_clt(filda_dict['FP_Latitude'].values,
										  filda_dict['FP_Longitude'].values,
										  nl, tt, **kwargs)
	for key in emis_dict.keys():
		filda_dict[key] = emis_dict[key]
	
	if nl.flag_emit:
		message = f'Setup to use NRT surface emissivity...'
		printf(message, 1, prefix = 'MCBEF_IO')
		emis_dict = get_surface_emit_sinu(filda_dict['FP_Latitude'].values,
										  filda_dict['FP_Longitude'].values,
										  nl, tt, **kwargs)
		for key in emis_dict.keys():
			filda_dict[key] = emis_dict[key]
		
		# Fill NaN values in 'Emis_xxx' with values from 'Emis_xxx_CLT'
		for param in emit_params:
			filda_dict[param] = filda_dict[param].fillna(filda_dict[param + '_CLT'])
	
	else:
		for param in emit_params:
			filda_dict[param] = filda_dict[param + '_CLT']
		
		
	return filda_dict


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_surface_emit_sinu(lat, lon, nl, tt, **kwargs):
 
	numCeil	  = kwargs.get('numCeil', 1200)
	
	effective_doy = str(int(float(tt.DOY)//8 * 8 + 1))
	
	DATA_DIR = nl.input_path + 'VNP21A2' + '/' +  tt.Y + '/' + effective_doy.zfill(3) + '/'
	
	tiles, hidMax, hidMin, vidMax, vidMin = get_tile_sinusoidal_3((lat, lon), numCeil = numCeil)
	
	GridDim = ((vidMax - vidMin + 1) * numCeil, (hidMax - hidMin + 1) * numCeil)
	
	params = ['Emis_14', 'Emis_15', 'Emis_16']
	data_grid_dict = {} 
	for param in params:
		data_grid_dict[param] = np.full(GridDim, np.nan)
	
	for tile in tiles:
		hid = int(float(tile[1:3]))
		vid = int(float(tile[4:]))
	
		hIdx = hid - hidMin 
		vIdx = vid - vidMin
	
		filename = glob.glob(DATA_DIR + '*' + tile + '*.h5')    
		if len(filename) == 0:
			continue
				
		if nl.flag_verbose:
			message = f'Reading {filename[0]}'
			printf(message, 1, prefix = 'MCBEF_IO')
	
		VNP21A2 = read_VNP21A2(filename[0], params)
	
		for param in params:
			data_grid_dict[param][vIdx * numCeil : (vIdx + 1) * numCeil, \
								  hIdx * numCeil : (hIdx + 1) * numCeil, ] = VNP21A2[param]
	
	
	boundary_tile = 'h' + str(hidMin).zfill(2) + 'v' + str(vidMin).zfill(2)
	x, y, resol = cal_sinu_xy(boundary_tile, numCeil)
	x_min = np.nanmin(x)
	y_max= np.nanmax(y)
	
	print(f' - Given {boundary_tile} X min: {x_min}, Y max: {y_max}, resolution: {resol}')
	
	x_s, y_s = geog_to_sinu([lat, lon])
	x_id = (x_s - x_min + resol/2.)//resol
	y_id = (y_max - y_s + resol/2.)//resol
	
	
	x_id = x_id.astype(int)
	y_id = y_id.astype(int)
	emis_dict = {}
	for param in params:
		band_name = param.split('_')[0] + '_M' + param.split('_')[1]
		emis_dict[band_name] = data_grid_dict[param][y_id, x_id]
	
	emis_dict['Emis_I05'] = 0.5*emis_dict['Emis_M15'] + 0.5 * emis_dict['Emis_M16']
		
	return emis_dict


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_surface_emit_sinu_clt(lat, lon, nl, tt, **kwargs):

	numCeil = kwargs.get('numCeil', 1200)
	
	effective_doy = str(int(float(tt.DOY)//8 * 8 + 1))
	
	# use the climatology path...
	DATA_DIR = nl.input_path + 'VNP21A2' + '/CLT/' + effective_doy.zfill(3) + '/'
	
	tiles, hidMax, hidMin, vidMax, vidMin = get_tile_sinusoidal_3((lat, lon), numCeil = numCeil)
	
	GridDim = ((vidMax - vidMin + 1) * numCeil, (hidMax - hidMin + 1) * numCeil)
	
	params = ['Emis_14', 'Emis_15', 'Emis_16', 'LST_Day_1KM', 
			  'LST_Night_1KM', 'LST_Day_1KM_STD', 'LST_Night_1KM_STD']
	
	data_grid_dict = {} 
	for param in params:
		data_grid_dict[param] = np.full(GridDim, np.nan)
	
	
	for tile in tiles:
		hid = int(float(tile[1:3]))
		vid = int(float(tile[4:]))
	
		hIdx = hid - hidMin 
		vIdx = vid - vidMin
		
		filename = glob.glob(DATA_DIR + '*' + tile + '*.nc')    
		if len(filename) == 0:
			continue
				
		if nl.flag_verbose:
			message = f'Reading {filename[0]}'
			printf(message, 1, prefix = 'MCBEF_IO')
	
		VNP21A2 = read_VNP21A2_CLT(filename[0], params)
	
		for param in params:
			data_grid_dict[param][vIdx * numCeil : (vIdx + 1) * numCeil, \
								  hIdx * numCeil : (hIdx + 1) * numCeil, ] = VNP21A2[param]    
	
	boundary_tile = 'h' + str(hidMin).zfill(2) + 'v' + str(vidMin).zfill(2)
	x, y, resol = cal_sinu_xy(boundary_tile, numCeil)
	x_min = np.nanmin(x)
	y_max= np.nanmax(y)
	
	x_s, y_s = geog_to_sinu([lat, lon])
	x_id = (x_s - x_min + resol/2.)//resol
	y_id = (y_max - y_s + resol/2.)//resol
	
	
	x_id = x_id.astype(int)
	y_id = y_id.astype(int)
	emis_dict = {}
	for param in params:
		if 'Emis' in param:
			band_name = param.split('_')[0] + '_M' + param.split('_')[1] + '_CLT'
			emis_dict[band_name] = data_grid_dict[param][y_id, x_id]
		else:
			emis_dict[param] = data_grid_dict[param][y_id, x_id]
	
	emis_dict['Emis_I05_CLT'] = 0.5*emis_dict['Emis_M15_CLT'] + 0.5 * emis_dict['Emis_M16_CLT']
	
	return emis_dict


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_VNP21A2(filename, params):
	import h5py
	import numpy as np
	output = {}
	datafiled = 'HDFEOS/GRIDS/VIIRS_Grid_8Day_1km_LST21/Data Fields/'
	print(' - Reading',filename)
	f = h5py.File(filename, 'r')
	
	# 	NO_SCLAE_PARAMS = ['Count_Day', 'Count_Night', 'QC_Day', 'QC_Night']
	
	for param in params:
		DNs = f[datafiled + param][:]
	
		scale_factor = f[datafiled + param].attrs['scale_factor'][0]
		add_offset   = f[datafiled + param].attrs['add_offset'][0]
		units        = f[datafiled + param].attrs['units']
		long_name    = f[datafiled + param].attrs['long_name'].decode('UTF-8')
	
		output[param] = DNs * 1.0 * scale_factor + add_offset
	
	
		valid_range  = f[datafiled + param].attrs['valid_range']
		vmin = valid_range[0]
		vmax = valid_range[1]		
		invalid = np.where((DNs<vmin) | (DNs>vmax))
	
		output[param][invalid] = np.nan
	
	f.close()
	
	return output

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def read_VNP21A2_CLT(filename, params):
    
    print(f' - Reading {filename}')
    output = {}
    ncid = Dataset(filename, 'r')
    ncid.set_auto_mask(False)
    ncid.set_auto_scale(False)
    
    for param in params:

        scale_factor = ncid[param].scale_factor
        add_offset = ncid[param].add_offset
        valid_range = ncid[param].valid_range
        fill_value  = ncid[param]._Fillvalue
        DNs = ncid[param][:].T
        invalid = np.where( (DNs < valid_range[0]) | \
                            (DNs > valid_range[1]) )

        output[param] = DNs.astype(float) * scale_factor + add_offset
        output[param][invalid] = np.nan

    ncid.close()

    return output

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_static_thermal_anomaly(filda_dict, nl, tt):
	'''
	'''
	DATA_DIR = nl.input_path + '/STATIC/static_thermal_anomaly.' + str(tt.Y) + '.nc'
	ncid = Dataset(DATA_DIR, 'r')
	ncid.set_auto_mask(False)
	static_mask = ncid['static_mask'][:]
	
	W = ncid.W
	S = ncid.S
	resol_lat = ncid.resol_lat
	resol_lon = ncid.resol_lon
	
	ncid.close()
	
	idx_lon = np.round((filda_dict['FP_Longitude'] - W )//resol_lon).astype(int)
	idx_lat = np.round((filda_dict['FP_Latitude'] - S)//resol_lat).astype(int)    
	
	filda_dict['Static_flag'] = static_mask[idx_lon, idx_lat]
	
	return filda_dict


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def read_MCD12Q1(filename, params, verbose = False):
	'''
	Read the MCD12Q1 with pyhdf API

	'''
	
	from pyhdf.SD import SD, SDC
	import numpy as np
	
	if verbose:
		print(f' - Reading', filename)
	
	# open the HDF4 file
	hdf = SD(filename, SDC.READ)
	
	output = {}
	
	for param in params:
		
		# get the Land_Cover_Type_1 data
		sds = hdf.select(param)
		data = sds.get()
		# get the valid range and fill value metadata
		valid_range = sds.attributes()['valid_range']
		fill_value = sds.attributes()['_FillValue']

		# apply the valid range mask and fill value mask
		data = np.ma.masked_outside(data, valid_range[0], valid_range[1])
		data = np.ma.masked_equal(data, fill_value)	

		output[param] = data
		
	return output

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_surface_type(filda_dict, nl, tt, defalt_year = None, **kwargs):

	if defalt_year is not None:
		year = str(int(defalt_year))
		message = f'Use the default land surface database [MODIS {year}].'
		printf(message, 1, prefix = 'MCBEF_IO')
	else:
		year = tt.Y
	
	DATA_DIR = nl.input_path + 'MCD12Q1' + '/' +  year + '/001/'

	land_type = get_surface_type_sinu(filda_dict['FP_Latitude'], filda_dict['FP_Longitude'], DATA_DIR)
	filda_dict['FP_Land_Type'] = land_type.astype(int)

	return filda_dict

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_surface_type_sinu(lat, lon, DATA_DIR, **kwargs):
	'''
	provide a set of coordinates, return the land surface type
	
	'''
	numCeil	  = kwargs.get('numCeil', 2400)

	tiles, hidMax, hidMin, vidMax, vidMin = get_tile_sinusoidal_3((lat, lon), numCeil = numCeil)
	
	GridDim = ((vidMax - vidMin + 1) * numCeil, (hidMax - hidMin + 1) * numCeil)
	
	params = ['LC_Type1']
	
	data_grid_dict = {} 
	for param in params:
		data_grid_dict[param] = np.full(GridDim, np.nan)
		
	for tile in tiles:
		hid = int(float(tile[1:3]))
		vid = int(float(tile[4:]))
	
		hIdx = hid - hidMin 
		vIdx = vid - vidMin
		filename = glob.glob(DATA_DIR + '*' + tile + '*.hdf')    
		if len(filename) == 0:
			continue
			
		MCD12Q1 = read_MCD12Q1(filename[0], params)

		for param in params:
			data_grid_dict[param][vIdx * numCeil : (vIdx + 1) * numCeil, \
								  hIdx * numCeil : (hIdx + 1) * numCeil, ] = MCD12Q1[param]

	boundary_tile = 'h' + str(hidMin).zfill(2) + 'v' + str(vidMin).zfill(2)
	x, y, resol = cal_sinu_xy(boundary_tile, numCeil)
	x_min = np.nanmin(x)
	y_max= np.nanmax(y)
	
	print(f' - Given {boundary_tile} X min: {x_min}, Y max: {y_max}, resolution: {resol}')
	x_s, y_s = geog_to_sinu([lat, lon])
	x_id = (x_s - x_min + resol/2.)//resol
	y_id = (y_max - y_s + resol/2.)//resol
	
	x_id = x_id.astype(int)
	y_id = y_id.astype(int)
	land_type = data_grid_dict['LC_Type1'][y_id, x_id]

	return land_type

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Function to write
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def ini_output_dir(nl, tt):

    '''
    ini_output_dir creates the out put path for a given time

    '''
    #-- it is no necessary to send the manelist dictionary, just send the out_dir path

    sample_path = nl.output_path + 'SAMPLE/'
    if os.path.isdir(sample_path) == False:
        message = f'Making directory {sample_path}'
        print(message)
        os.mkdir(sample_path)

    sample_path = sample_path + tt.Y + '/'
    if os.path.isdir(sample_path) == False:
        message = f'Making directory {sample_path}'
        print(message)
        os.mkdir(sample_path)

    sample_path = sample_path + tt.DOY + '/'
    if os.path.isdir(sample_path) == False:
        message = f'Making directory {sample_path}'
        print(message)
        os.mkdir(sample_path)


    state_path = nl.output_path + 'STATE/'
    if os.path.isdir(state_path) == False:
        message = f'Making directory {state_path}'
        print(message)
        os.mkdir(state_path)

    state_path = state_path + tt.Y + '/'
    if os.path.isdir(state_path) == False:
        message = f'Making directory {state_path}'
        print(message)
        os.mkdir(state_path)

    state_path = state_path + tt.DOY + '/'
    if os.path.isdir(state_path) == False:
        message = f'Making directory {state_path}'
        print(message)
        os.mkdir(state_path)


    return sample_path, state_path

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def dict2csv(dataDic, filename, badvalue = 'N/A', verbose=0):
	'''
	'''
	
	message = f'Writing {filename}...'
	printf(message, verbose, prefix = 'MCBEF_IO')
	
	df = pd.DataFrame.from_dict(dataDic,'columns')
	df.to_csv(filename, index=False, na_rep = 'N/A')
	return df

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def parse_keywords_output(file_path):
	keywords = {}
	with open(file_path, 'r') as file:
		for line in file:
			 # Skip empty lines and comments
			if line.strip() and not line.startswith('#'):
				variable, long_name, unit, legend, data_type  = line.strip().split(';')
				keywords[variable.strip()] = [long_name.strip(), unit.strip(), legend.strip(), data_type.strip()]
				
	return keywords

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def dict2nc(save_dict, nl, sat, data_type, savename):

	'''
	write_nc writes the provided dictionary into *.nc files
	
	'''
	
	# Get the absolute path to the directory where the current 
	base_dir = os.path.dirname(os.path.abspath(__file__))
	# Append the relative path to the desired directory
	output_dir = os.path.join(base_dir, 'output.rc')
	# Defined at a module or global level	
	infor_dict = parse_keywords_output(output_dir)

	# in case some of the parameters are not in numpy array 
	for key in save_dict.keys():
		save_dict[key] = np.array(save_dict[key])

	if data_type == 'Sample':
		fire, sample = np.shape(save_dict['t_s'])
	
		if fire <=0:
			print(' - MCBEF: No fire detected...')
			return 	
		ncid = Dataset(savename, 'w', format='NETCDF4' )
		ncid.createDimension('fire', fire)
		ncid.createDimension('sample', sample)
	
	if data_type == 'State':
		fire = save_dict['t_s_upp'].shape[0]
		if fire <=0:
			print(' - MCBEF: No fire detected...')
			return 	
		ncid = Dataset(savename, 'w', format='NETCDF4' )
		ncid.createDimension('fire', fire)
	
	# Define a CRS variable for the geographic coordinate system
	crs_var = ncid.createVariable('crs', 'i4')
	crs_var.grid_mapping_name = "latitude_longitude"
	crs_var.epsg_code = "EPSG:4326"  # WGS84 standard
	crs_var.semi_major_axis = 6378137.0
	crs_var.inverse_flattening = 298.257223563
	crs_var.longitude_of_prime_meridian = 0.0
		
		
	for key in save_dict.keys():
	
		dataType = type(save_dict[key])
		if dataType == np.ndarray:
			nDim = len(np.shape(save_dict[key]))
			if nDim == 1:
				if key in infor_dict.keys():
		 
					dtype = infor_dict[key][3]
					tempInstance = ncid.createVariable( key, dtype, ('fire'), zlib=True, complevel = 4 )
					tempInstance[:] = save_dict[key][:]
	
					tempInstance.units = infor_dict[key][1]
					tempInstance.long_name = infor_dict[key][0]
					tempInstance.legend    = infor_dict[key][2]
					tempInstance.data_type = dtype
					if key == 'FP_Longitude':
						tempInstance.standard_name = 'longitude'
					elif key == 'FP_Latitude':
						tempInstance.standard_name = 'latitude'
					else:
						# not necessary for Geographic Coordinates (Lat/Lon), but critial for
						# Sinusodial, UTM or other complex projections
						# keep it here for climate and forecast (CF) conventions
						tempInstance.grid_mapping = 'crs'
						tempInstance.coordinates = 'FP_Longitude FP_Latitude'
						
			if nDim == 2:
				if key in infor_dict.keys():
					dtype = infor_dict[key][3]
					tempInstance = ncid.createVariable( key, dtype, ('fire', 'sample'), zlib=True, complevel = 9) #, shuffle=True, chunksizes = (32, 6400)
					tempInstance[:, :] = save_dict[key][:, :]
	
					tempInstance.units     = infor_dict[key][1]
					tempInstance.long_name = infor_dict[key][0]
					tempInstance.legend    = infor_dict[key][2]
					tempInstance.data_type = dtype
	
	metadata_dict={}
	metadata_dict['satellite'] = sat
	
	for key in nl.__dict__.keys():
		metadata_dict[key] = getattr(nl, key)                    
					
	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	ncid.LongName = 'MCBEF fire parameter estimation on FILDA ' + metadata_dict['satellite'] + '47MOD.'
	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	
	for key, value in metadata_dict.items():
		# Convert booleans to integers (1 for True, 0 for False)
		if isinstance(value, bool):
			value = int(value)
		# Convert lists and dictionaries to strings because NetCDF attributes need to be strings or numbers
		elif isinstance(value, (list, dict)):
			value = str(value)
		# Set each key-value pair as a global attribute
		setattr(ncid, key, value)    
	
	ncid.close()
	print(' - MCBEF: Saving ' + savename)
	
	return



