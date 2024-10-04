'''
Comment in progress
'''
import FILDA_Time_Cord
import FILDA_IO


import os
import numpy as np
import glob
#-----------------------------------------------------------------------
def get_dnb_clt(modData, clt_var, time, namelist, verbose = False): 
	'''
	Function to prepare the DNB climatologies...
	-----------
	Parameters:
	
	latitude  : array, latitude of the region
	longitude : array, longitude of the region
	clt_var   : list, contains name of the climatologies to be generated
	year	  : string, years of the of climatologies 	
	month	  : string, years of the of climatologies 	
	
	-----------
	Returns: 
	output	  : dict, contains the climatologies
	
	'''

	
	### import os
	### import numpy as np

	# MZ, the prefix of the peatland now can be specified by the namelist...
	NTL_PREFIX = namelist['NTL_prefix']
	if NTL_PREFIX == 'N/A':
		print(f' - get_peatland: No prefix for peatland is specified...')
		NTL_PREFIX = ''	
	
	print(' - FILDA: Get DNB climatologies')
	# MZ change all the directory as the global variable
	NTL_DIR= namelist['NTL_DIR']
	north = np.nanmax(modData['latitude'])
	south = np.nanmin(modData['latitude'])
	west  = np.nanmin(modData['longitude'])
	east  = np.nanmax(modData['longitude'])
	cord  = [north, south, west, east]

	numCeil 	= 2400
	resol 		= 10./numCeil
	resol_half	= 10./numCeil/2.
	
	tiles = FILDA_Time_Cord.get_tiles(cord)
# 	print(f' - FILDA2: Need tiles:{tiles}')
	meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = FILDA_Time_Cord.get_cord_PlateCarree(cord, numCeil)
	lat_max = np.max(meshLat)
	lat_min = np.min(meshLat)
	
	lon_max = np.max(meshLon)
	lon_min = np.min(meshLon)

	latIdx = ( lat_max  - modData['latitude']  )//resol
	lonIdx = ( modData['longitude'] - lon_min )//resol
	
	invalid_idx = np.where( (modData['latitude'] != modData['latitude']) | (modData['longitude'] != modData['longitude']) )
	latIdx[invalid_idx] = 0
	lonIdx[invalid_idx] = 0	
	
	latIdx = latIdx.astype(int)
	lonIdx = lonIdx.astype(int)
	
	# create the climatology variables...
	output = {}
	for var in clt_var:
		output[var] = np.full_like(meshLat, np.nan)

	if float(time.Y) >= namelist['Default_year_NTL']:
		NTL_year = str(int(namelist['Default_year_NTL']))
	else:
		NTL_year = time.Y

	
# 	NTL_JDN = FILDA_Time_Cord.JulianDay(year+time.M+'01', outtype = 'nasa')
	
	NTL_JDN = FILDA_Time_Cord.JulianDay(NTL_year+time.M+'01', outtype = 'nasa')
	
	i = 0
	for tile in tiles:
		hh = np.int( np.float( tile[1:3] ) )
		vv = np.int( np.float( tile[4:]  ) )
		
		# MZ: switch to a flexible way to find filename
		# MZ May-11-2023, modify to standard name for the NTL
# 		NTL_DIR = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/NIGHTTIME_LIGHT/THREE_MONTH_AVG/'
# 		year = '2019'
# 		strMonth = '09'
# 		filename = [NTL_DIR + year + strMonth + '.' + tile + '.nc']
# 		if os.path.exists(filename[0]):
		filename = glob.glob( NTL_DIR + '*' + NTL_PREFIX + '*' + NTL_JDN + '.' + tile + '*.nc' )
		if len(filename)>0:
			nc_data    =  FILDA_IO.read_nc(filename[0], clt_var, verbose)
			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			for key in nc_data.keys():
				output[key][vIdx * numCeil : (vIdx + 1) * numCeil, \
					 		hIdx * numCeil : (hIdx + 1) * numCeil, ] = nc_data[key]
			i = i + 1
		else:
			pass
# 			print(' - Warning! get_dnb_clt: CANNOT FIND', NTL_JDN + '.' + tile + '.nc for lighttime light')

	for key in output.keys():
		output[key]  = output[key][latIdx, lonIdx]
		output[key][latIdx[invalid_idx],lonIdx[invalid_idx]] = np.nan
	
	return output
#-----------------------------------------------------------------------
def get_dnb_clt_2(latitude, longitude, clt_var, year, month):
	'''
	Function to prepare the DNB climatologies, currently 4 month of the 
	nighttime light is used to derive the nighttime light database...
	
	-----------
	Parameters:
	
	latitude  : array, latitude of the region
	longitude : array, longitude of the region
	clt_var   : list, contains name of the climatologies to be generated
	year	  : string, years of the of climatologies 	
	month	  : string, years of the of climatologies 	
	
	-----------
	Returns: 
	output	  : dict, contains the climatologies
	
	'''
	### import os
	### import numpy as np
	
	print(' - get DNB climatologies')
	
	NTL_DIR= '/Dedicated/jwang-data/shared_satData/VIIRS/VNP46A1/climatologies/tile/'
	north = np.max(latitude)
	south = np.min(latitude)
	west  = np.min(longitude)
	east  = np.max(longitude)
	cord  = [north, south, west, east]
	
	numCeil = 2400
	resol = 10./numCeil
	tiles = FILDA_Time_Cord.get_tiles(cord)
	meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = FILDA_Time_Cord.get_cord_PlateCarree(cord, numCeil)
	
	lat_max = np.max(meshLat)
	lat_min = np.min(meshLat)
	
	
	lon_max = np.max(meshLon)
	lon_min = np.min(meshLon)
	
	
	latIdx = (lat_max - latitude)//resol
	latIdx = latIdx.astype(int)
	lonIdx = (longitude -lon_min)//resol
	lonIdx = lonIdx.astype(int)
	

	output = {}
		
	for var in clt_var:

		std_var = np.full_like(meshLat, np.nan)
		mean_var = np.full_like(meshLat, np.nan)
		
		varnames = []
		month = int( float(month) )
		month = np.arange( month-5, month, 1)
		months = []
		for m in month:
			if m<1:
				m = m + 12
			if m == 1:
				m = 2
			strMonth = str(m)
			if len(strMonth)<2:
				strMonth = '0' + strMonth	
			months.append(strMonth)

		varnames = [var+'_'+year+month for month in months]
		
		
		for tile in tiles:
			hh = np.int( np.float( tile[1:3] ) )
			vv = np.int( np.float( tile[4:]  ) )

			filename_mean = NTL_DIR + 'VNP46A1.' + tile + '.' + var + '.mean.nc'
			filename_std  = NTL_DIR + 'VNP46A1.' + tile + '.' + var + '.std.nc'
	
			if os.path.exists(filename_mean):

				mean_data    =  FILDA_IO.read_nc( filename_mean, varnames)
				std_data     =  FILDA_IO.read_nc( filename_std,  varnames)
							
				mean_data_array = []
				for key in mean_data:
					mean_data_array.append(mean_data[key])
				
				std_data_array = []
				for key in std_data:
					std_data_array.append(std_data[key])
					
				mean_data_array = np.stack(mean_data_array, axis = 2)
				mean_data_array = np.nanmean(mean_data_array, axis = 2)

				std_data_array = np.stack(std_data_array, axis = 2)
				std_data_array = np.nanmean(std_data_array, axis = 2)
				
				
				hIdx = hh - hidMin 
				vIdx = vv - vidMin
				mean_var[vIdx * numCeil : (vIdx + 1) * numCeil, \
						 hIdx * numCeil : (hIdx + 1) * numCeil, ] = mean_data_array
						
				std_var[vIdx * numCeil : (vIdx + 1) * numCeil, \
						hIdx * numCeil : (hIdx + 1) * numCeil, ]  = std_data_array

		output[var + '_std']  = std_var[latIdx, lonIdx]
		output[var + '_mean'] = mean_var[latIdx, lonIdx]

	return output
