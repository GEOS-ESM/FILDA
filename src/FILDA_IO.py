'''

Python library of FILDA fire detection


IO module

This module provide the basice reading and writing functions for FILDA

1. READ_IMG

2. READ_MOD

3. READ_DNB

4. read_viirs_i_band 

5. read_viirs_m_band

6. read_viirs_dnb

7. read_Geos_FP

8. write_nc_img

9. write_nc_mod

# update in Feb.28 2022, use the build in BT lut to calculate BT or each channel

'''
import numpy as np
import copy
import glob
from netCDF4 import Dataset
import sys, os
import time
import pandas as pd
import FILDA_Time_Cord
import shutil
#-----------------------------------------------------------------------
def get_files(sat, inputDir, time, mode = 1):
	
	if mode == 1:
		file_dict = get_files_mod_1(sat, inputDir, time)
	if mode == 2:
		file_dict = get_files_mod_2(sat, inputDir, time)
	if mode == 3:
		file_dict = get_files_mod_3(sat, inputDir, time)
	print(f' - FILDA_IO: Reading data through mode {mode}')
	return file_dict


#-----------------------------------------------------------------------
def get_files_mod_3(sat, namelist, time):
	'''
	
	Function to pair the level-B and the geos-fp data together Mode 3.
	
	This mode match the level 1b file and geos-fp files in a redefined 
	RUN directory.


	'''

	# put all files into a dictionary
	file_dict = {}
	found_files = True
	#-------------
	# DNB radiance
	#-------------

	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '02DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	try:
		file_dict[sat + '02DNB'] = file_temp[0]
		DayNightFlag  = read_DayNightFlag(file_temp[0])
		if DayNightFlag == 'Day':
			file_dict['DayNightFlag'] = 'Day'
			return file_dict
		file_dict['DayNightFlag'] = 'Night/Both' 
	
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '02DNB')
		found_files = False
	

	#----------------
	# DNB geolocation
	#----------------
	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '03DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	
	try:
		file_dict[sat + '03DNB'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '03DNB')
		found_files = False

	#-------------
	# MOD radiance
	#-------------	
	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '02MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')
	
	try:
		file_dict[sat + '02MOD'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '02MOD')
		found_files = False
		
	#----------------
	# MOD geolocation
	#----------------	
	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '03MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '03MOD'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '03MOD')
		found_files = False

		
	#-------------
	# IMG radiance
	#-------------			
	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '02IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '02IMG'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '02IMG')
		found_files = False

	#----------------
	# IMG geolocation
	#----------------				
	file_temp = glob.glob( namelist['RUN_DIR'] + sat + '03IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '03IMG'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find", sat + '03IMG')
		found_files = False
	

	#-------------
	# GEOS-FP
	#-------------
	file_temp = glob.glob( namelist['RUN_DIR'] + '*' + time.GOES + '*')

	try:
		file_dict['GEOS-FP'] = file_temp[0]
	except (IndexError):
		print(" - FILDA_IO: Could not find GEOS-FP")
		found_files = False
	
	if not found_files:
		file_dict = {}
	

	return file_dict


#-----------------------------------------------------------------------
def get_files_mod_1(sat, namelist, time):

	'''
	Function to pair the level-B and the geos-fp data together Mode 1.
	
	This mode match the level 1b file and geos-fp files that are stored in 
	NASA LADDS layout.

	'''
	### import glob
	# Directory of DNB radiance data 
	dnbRadDir  = namelist['RUN_DIR'] + sat + '02DNB/' + time.Y + '/' + time.DOY + '/'
	# Directory of DNB geo data
	dnbGeoDir  = namelist['RUN_DIR'] + sat + '03DNB/' + time.Y + '/' + time.DOY + '/'
	# Directory of Mband radiance data
	modRadDir  = namelist['RUN_DIR'] + sat + '02MOD/' + time.Y + '/' + time.DOY + '/'
	# Directory of Mband geo data
	modGeoDir  = namelist['RUN_DIR'] + sat + '03MOD/' + time.Y + '/' + time.DOY + '/'
	# Directory of I band['INPUT_DIR'] radiance data
	imgRadDir  = namelist['RUN_DIR'] + sat + '02IMG/' + time.Y + '/' + time.DOY + '/'
	# Directory of i band geo data
	imgGeoDir  = namelist['RUN_DIR'] + sat + '03IMG/' + time.Y + '/' + time.DOY + '/'
	# Directory of geos-fp
	geosFPDir  = namelist['RUN_DIR'] + 'GEOS_FP/'

	# put all files into a dictionary
	file_dict = {}
	found_files = True
	#-------------
	# DNB radiance
	#-------------
	
	file_temp = glob.glob( dnbRadDir + sat + '02DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	
	try:
		file_dict[sat + '02DNB'] = file_temp[0]
		DayNightFlag  = read_DayNightFlag(file_temp[0])
		if DayNightFlag == 'Day':
			file_dict['DayNightFlag'] = 'Day'
			return file_dict
		file_dict['DayNightFlag'] = 'Night/Both' 
	
	except (IndexError):
		print("Could not find", sat + '02DNB')
		found_files = False
	
	#----------------
	# DNB geolocation
	#----------------
	file_temp = glob.glob( dnbGeoDir + sat + '03DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	
	try:
		file_dict[sat + '03DNB'] = file_temp[0]
	except (IndexError):
		print("Could not find", sat + '03DNB')
		found_files = False

	#-------------
	# MOD radiance
	#-------------	
	file_temp = glob.glob( modRadDir + sat + '02MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')
	
	try:
		file_dict[sat + '02MOD'] = file_temp[0]
	except (IndexError):
		print("Could not find", sat + '02MOD')
		found_files = False
		
	#----------------
	# MOD geolocation
	#----------------	
	file_temp = glob.glob( modGeoDir + sat + '03MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '03MOD'] = file_temp[0]
	except (IndexError):
		print("Could not find", sat + '03MOD')
		found_files = False

		
	#-------------
	# IMG radiance
	#-------------			
	file_temp = glob.glob( imgRadDir + sat + '02IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '02IMG'] = file_temp[0]
	except (IndexError):
		print("Could not find", sat + '02IMG')
		found_files = False

	#----------------
	# IMG geolocation
	#----------------				
	file_temp = glob.glob( imgGeoDir + sat + '03IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')

	try:
		file_dict[sat + '03IMG'] = file_temp[0]
	except (IndexError):
		print("Could not find", sat + '03IMG')
		found_files = False
	

	#-------------
	# GEOS-FP
	#-------------
	file_temp = glob.glob( geosFPDir + '*' + time.GOES + '*')

	try:
		file_dict['GEOS-FP'] = file_temp[0]
	except (IndexError):
		print("Could not find GEOS-FP")
		found_files = False
	
	if not found_files:
		file_dict = {}
	
	return file_dict


#-----------------------------------------------------------------------
def get_files_mod_2(sat, namelist, time):

	### import glob
	# Directory of DNB radiance data 
	dnbRadDir  = namelist['INPUT_DIR'] + sat + '02DNB_N/'
	# Directory of DNB geo data
	dnbGeoDir  = namelist['INPUT_DIR'] + sat + '03DNB_N/'
	# Directory of Mband radiance data
	modRadDir  = namelist['INPUT_DIR'] + sat + '02MOD_N/'
	# Directory of Mband geo data
	modGeoDir  = namelist['INPUT_DIR'] + sat + '03MOD_N/'
	# Directory of I band radiance data
	imgRadDir  = namelist['INPUT_DIR'] + sat + '02IMG_N/'
	# Directory of i band geo data
	imgGeoDir  = namelist['INPUT_DIR'] + sat + '03IMG_N/'
	# Directory of geos-fp
	geosFPDir  = namelist['INPUT_DIR'] + 'geos_fp/'

	# put all files into a dictionary
	file_dict = {}
	#-------------
	# DNB radiance
	#-------------
	file_temp = glob.glob( dnbRadDir + sat + '02DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '02DNB')
	else:
		file_dict[sat + '02DNB'] = file_temp[0]
	
	
	DayNightFlag  = read_DayNightFlag(file_temp[0])
	if DayNightFlag == 'Day':
		file_dict['DayNightFlag'] = 'Day'
		return file_dict
	file_dict['DayNightFlag'] = 'Night/Both'
	
	#----------------
	# DNB geolocation
	#----------------
	file_temp = glob.glob( dnbGeoDir + sat + '03DNB' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '03DNB')
	else:
		file_dict[sat + '03DNB'] = file_temp[0]

	#-------------
	# MOD radiance
	#-------------	
	file_temp = glob.glob( modRadDir + sat + '02MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '02MOD')
	else:
		file_dict[sat + '02MOD'] = file_temp[0]
		
	#----------------
	# MOD geolocation
	#----------------	
	file_temp = glob.glob( modGeoDir + sat + '03MOD' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '03MOD')
	else:
		file_dict[sat + '03MOD'] = file_temp[0]
		
	#-------------
	# IMG radiance
	#-------------			
	file_temp = glob.glob( imgRadDir + sat + '02IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '02IMG')
	else:
		file_dict[sat + '02IMG'] = file_temp[0]

	#----------------
	# IMG geolocation
	#----------------				
	file_temp = glob.glob( imgGeoDir + sat + '03IMG' + '.' + time.JDN + '.' + time.OP + '*.nc')
	if len(file_temp)<= 0:
		raise Exception("Could not find", sat + '03IMG')
	else:
		file_dict[sat + '03IMG'] = file_temp[0]

	#-------------
	# GEOS-FP
	#-------------
	file_temp = glob.glob( geosFPDir + '*' + time.GOES + '*')
	if len(file_temp)<= 0:
		print(geosFPDir + '*' + time.GOES + '*')
		raise Exception("Could not find GEOS-FP")
	else:
		file_dict['GEOS-FP'] = file_temp[0]

	
	
	
	return file_dict


#-----------------------------------------------------------------------
def READ_IMG(imgRadFile, imgGeoFile, namelist, remove_flag = True):
	'''
	High level function for reading the VIIRS IMG data...
	'''
	### import numpy as np
	### import copy
	#-------------------------------
	#read I band data...
	#-------------------------------
	imgVars = ['I04', 'I05', 'I04_quality_flags', 'I05_quality_flags']
# 	imgGeoVars = ['latitude', 'longitude', 'land_water_mask','sensor_zenith']
	imgGeoVars = ['latitude', 'longitude', 'land_water_mask', 'solar_zenith', 
	              'sensor_zenith', 'quality_flag', 'sensor_azimuth']
	imgData = read_viirs_i_band(radFile = imgRadFile, 
								radVars = imgVars,
								geoFile = imgGeoFile,
								geoVars = imgGeoVars, )
	

	#-------------------------------							
	# filter by solar zenith and latitude, polar region are exclude from testing	
	invalidIdx = np.where( (imgData['I04_quality_flags']>=128) |\
						   (imgData['I05_quality_flags']>=128) |\
						   (imgData['solar_zenith']<namelist['twilight_ang_min']) |\
						   (imgData['latitude']<-60)           |\
						   (imgData['latitude']>85)            |\
						   (imgData['quality_flag']>0)         )
	qaPara = ['I04_rad', 'I05_rad']
	for para in qaPara:
		imgData[para][invalidIdx] = np.nan


	imgData['land_water_mask_2'] = copy.deepcopy(imgData['land_water_mask'])
	land_water_mask_2 = np.ones_like(imgData['land_water_mask'])
	# add 6 = Continental for water masking
	idx               = np.where((imgData['land_water_mask']==0) | \
								 (imgData['land_water_mask']==7) | \
								 (imgData['land_water_mask']==6) )
# 	idx               = np.where((imgData['land_water_mask']!=1) | (imgData['land_water_mask']!=1) )
	land_water_mask_2[idx] = 0
	land_water_mask_2 = land_water_mask_2.astype(int)
	imgData['land_water_mask'] = land_water_mask_2
	
	# MZ 04-23-2024
	# enlarge SAA domain
	cord_SAA = [7, -55, -110, 20]
	
	idx = np.where( (imgData['latitude']  < cord_SAA[0])	& \
					(imgData['latitude']  > cord_SAA[1])	& \
					(imgData['longitude'] > cord_SAA[2])	& \
					(imgData['longitude'] < cord_SAA[3]))
					
	ssa_flag = np.zeros_like(imgData['latitude']).astype(int)
	ssa_flag[idx] = 1
	imgData['SAA_flag'] = ssa_flag
	
	return imgData
	
#-----------------------------------------------------------------------
def READ_MOD(modRadFile, modGeoFile, namelist, remove_flag = True):
	'''
	High level function for reading the VIIRS MOD data...	
	'''
	### import numpy as np
	### import copy
	#-------------------------------
	# read m band data...
	#-------------------------------
	# nighttime MOD does not have M09...
	
	modVars = ['M07', 'M08', 'M10', 'M11', 'M12', 'M13', 
	           'M14', 'M15', 'M16', 'M13_quality_flags', 
	           'M15_quality_flags','M11_quality_flags', 'M10_quality_flags']
	mGeoVars = ['latitude', 'longitude', 'land_water_mask','solar_zenith', 'sensor_zenith', 'sensor_azimuth']
	modData = read_viirs_m_band( radFile = modRadFile, 
							     radVars = modVars,
							     geoFile = modGeoFile,
							     geoVars = mGeoVars, 
							     navigation = False)
							     
	#-------------------------------							
	# filter out the invalid pixel
	#-------------------------------	
	# zero value does not physically make sense..
	invalidIdx_mod = np.where(  (modData['M11_quality_flags'] == 4) |  (modData['M11_quality_flags'] >= 128))
	qaPara = ['M07_rad', 'M08_rad', 'M10_rad', 'M11_rad']
	for para in qaPara:
		modData[para][invalidIdx_mod] = np.nan

	# first by bad radiance	
	qaPara = ['M07_rad', 'M08_rad', 'M10_rad', 'M11_rad']
	for para in qaPara:	
		invalidIdx_mod = np.where( (modData[para]<=0) )
		modData[para][invalidIdx_mod] = np.nan
	
	
	# second by QA flags
	invalidIdx_mod = np.where( (modData['M13_quality_flags'] == 4) |  (modData['M13_quality_flags'] >= 128) ) # | (modData['M15_quality_flags']>0)
	qaPara = ['M13_rad', 'M14_rad', 'M15_rad', 'M16_rad']
	for para in qaPara:
		modData[para][invalidIdx_mod] = np.nan
	
	# filter by solar zenith and latitude, polar region are exclude from testing
	invalidIdx_mod = np.where( (modData['solar_zenith']<namelist['twilight_ang_min']) | (modData['latitude']<-60) | (modData['latitude']>85))
	qaPara = ['M07_rad', 'M08_rad', 'M10_rad', 'M11_rad', 'M13_rad', 'M14_rad', 'M15_rad', 'M16_rad']
	for para in qaPara:
		modData[para][invalidIdx_mod] = np.nan
	
	invalidIdx_mod = np.where(modData['M13_rad']!=modData['M13_rad'])
	
	if remove_flag:
		removed_var = ['M11_quality_flags', 'M10_quality_flags',]
		for var in removed_var:
			del modData[var]
	# add 6 = Continental for water masking
	modData['land_water_mask_2'] = copy.deepcopy(modData['land_water_mask'])
	land_water_mask_2 = np.ones_like(modData['land_water_mask'])
	idx               = np.where((modData['land_water_mask']==0) | \
	                             (modData['land_water_mask']==7) | \
	                             (modData['land_water_mask']==6) )
# 	idx               = np.where((modData['land_water_mask']!=1) )
	land_water_mask_2[idx] = 0 
	land_water_mask_2 = land_water_mask_2.astype(int)
	modData['land_water_mask'] = land_water_mask_2
	
	# MZ 04-23-2024
	# enlarge SAA domain
	cord_SAA = [7, -55, -110, 20]
	
	idx = np.where( (modData['latitude']  < cord_SAA[0])	& \
					(modData['latitude']  > cord_SAA[1])	& \
					(modData['longitude'] > cord_SAA[2])	& \
					(modData['longitude'] < cord_SAA[3]))
					
	ssa_flag = np.zeros_like(modData['latitude']).astype(int)
	ssa_flag[idx] = 1
	modData['SAA_flag'] = ssa_flag

	
	return modData, invalidIdx_mod

#-----------------------------------------------------------------------
def READ_DNB(dnbRadFile, dnbGeoFile, namelist):
	'''
	High level function for reading the VIIRS DNB data...
	'''
	### import numpy as np
	### import copy
	#-------------------------------
	# read DNB band data...
	#-------------------------------
	dnbVar		= ['DNB_observations', 'DNB_quality_flags']
	dnbGeoVar	= ['latitude', 'longitude', 'solar_zenith']

	dnbData     = read_viirs_dnb(radFile = dnbRadFile, 
							     radVars = dnbVar,
							     geoFile = dnbGeoFile,
							     geoVars = dnbGeoVar, 
							     navigation = False)
	
	# filter by solar zenith and latitude, polar region are exclude from testing
	invalidIdx_dnb = np.where( (dnbData['solar_zenith']<namelist['twilight_ang_min']) | (dnbData['latitude']<-60) | (dnbData['latitude']>85))
	
	dnbData['DNB_observations'][invalidIdx_dnb] = np.nan
	# unit convert from Watts/cm^2/steradian tp nW/cm^2/steradian
	dnbData['DNB_observations'] = dnbData['DNB_observations']*1e9
	
	if 'VJ1' in dnbRadFile:
		for key in dnbVar + dnbGeoVar:
			dnbData[key] = dnbData[key][:, 17:3774]
			
	# MZ 04-23-2024
	# enlarge SAA domain
	cord_SAA = [7, -55, -110, 20]
	
	idx = np.where( (dnbData['latitude']  < cord_SAA[0])	& \
					(dnbData['latitude']  > cord_SAA[1])	& \
					(dnbData['longitude'] > cord_SAA[2])	& \
					(dnbData['longitude'] < cord_SAA[3]))
	
	ssa_flag = np.zeros_like(dnbData['latitude']).astype(int)
	ssa_flag[idx] = 1
	dnbData['SAA_flag'] = ssa_flag

	return dnbData

#-----------------------------------------------------------------------
def read_viirs_i_band(**kwargs):
	'''
	Function to read VIIRS M-band file
	
	'''
	### import numpy as np
	### from netCDF4 import Dataset
	
	print(' - Calling read_viirs_m_band...\n')
	# path of the rad file
	radFile	= kwargs.get('radFile', None)
	# list of rad variable to be read
	radVars	= kwargs.get('radVars', None)
	# path of the geo file
	geoFile	= kwargs.get('geoFile', None)
	# list of geo variables to be read
	geoVars	= kwargs.get('geoVars', None)
	# logical switch to determine whether return the navigation data
	navigation	= kwargs.get('navigation', True)
	
	non_refl_list = ['I04', 'I05']
	outDict = {}
	
	if radFile is not None:
		if len(radVars) > 0:
		
			print(' - Reading', radFile, '\n')
			
			ncid = Dataset(radFile, 'r')
			varList = list(ncid['/observation_data'].variables)

			# also read some global attributes of the rad-file.
			outDict['GRingPointLatitude'] = ncid.time_coverage_end
			outDict['GRingPointLongitude'] = ncid.time_coverage_start		
			
			ncid.set_auto_scale(False)			
			# netCDF will automatically apply the scale_factor and the add_offset
			# to the variable if its metadata contains these two.
			# However, it will also apply to the invalid value or scale the 
			# for some variable, the conversion equation is not the same
			# with the build in one...
# 			ncid.set_auto_mask(False)
		
			for key in radVars:
				if key in varList:
					if 'quality_flags' in key:
						outDict[key] = ncid["/observation_data/" + key][:]
					elif 'uncert_index' in key:
						ucidx = ncid["/observation_data/" + key][:].astype(float)
						scale_factor = ncid["/observation_data/" + key].scale_factor
						# valid_max    = ncid["/observation_data/" + key].valid_max
						# valid_min    = ncid["/observation_data/" + key].valid_min
						# use valid max and min to filter out invalid value
						# ucidx[np.where((ucidx > valid_max) | (ucidx < valid_min))] = np.nan		
						uc			 = 1 + scale_factor * ucidx**2
						uc			 = np.ma.filled(uc, fill_value = np.nan)	
						outDict[key] = uc
					else:
						if key in non_refl_list:
						
							DN           = ncid["/observation_data/" + key][:]			
							lut = ncid["/observation_data/" + key + '_brightness_temperature_lut'][:]
							scale_factor = ncid["/observation_data/" + key].scale_factor
							add_offset   = ncid["/observation_data/" + key].add_offset
							
							rad		 	 = DN.astype(float)  * scale_factor + add_offset
							rad		     = np.ma.filled(rad, fill_value = np.nan)
							

														
							outDict[key + '_rad'] = rad
							outDict[key + '_LUT'] = lut[:]
							outDict['BT'+key ]    = lut[DN]
							outDict['BT'+key ][np.where( outDict[key + '_rad'] != outDict[key + '_rad']) ] = np.nan
							outDict['BT'+key ]    = np.ma.filled(outDict['BT'+key ], fill_value = np.nan)
												
							
						if key not in non_refl_list:
							refl = ncid["/observation_data/" + key][:].astype(float)
							scale_factor = ncid["/observation_data/" + key].scale_factor
							add_offset   = ncid["/observation_data/" + key].add_offset
							refl		 = refl * scale_factor + add_offset
							refl		 = np.ma.filled(refl, fill_value = np.nan)
							outDict[key + '_refl'] = refl							
						
							rad = ncid["/observation_data/" + key][:].astype(float)
							scale_factor = ncid["/observation_data/" + key].radiance_scale_factor
							add_offset   = ncid["/observation_data/" + key].radiance_add_offset
							rad		 	 = rad + add_offset * scale_factor 
							rad		 	 = np.ma.filled(rad, fill_value = np.nan)
							outDict[key + '_rad'] = rad	
							
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')
			ncid.close()
	
	
	if geoFile is not None:
		if len(geoVars) > 0:
			ncid = Dataset(geoFile, 'r')
			geoList = list(ncid['/geolocation_data'].variables)
			print(' - Reading', geoFile, '\n')
			
			# also read some global attributes of the geo-file.
			outDict['GRingPointLatitude'] = ncid.GRingPointLatitude
			outDict['GRingPointLongitude'] = ncid.GRingPointLongitude
		
			outDict['SouthBoundingCoordinate'] = ncid.SouthBoundingCoordinate
			outDict['NorthBoundingCoordinate'] = ncid.NorthBoundingCoordinate
			outDict['WestBoundingCoordinate'] = ncid.WestBoundingCoordinate
			outDict['EastBoundingCoordinate'] = ncid.EastBoundingCoordinate

			outDict['DayNightFlag'] = ncid.DayNightFlag
			outDict['startDirection'] = ncid.startDirection
			outDict['endDirection'] = ncid.endDirection
			
			for key in geoVars:
				if key in geoList:
					temp = ncid["/geolocation_data/" + key][:]	
					outDict[key] = np.ma.filled(temp.astype(float), fill_value = np.nan)
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')

			if navigation:
				temp = ncid["/navigation_data/earth_moon_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_moon_distance'] = temp

				temp = ncid["/navigation_data/earth_sun_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_sun_distance'] = temp

			ncid.close()

	return outDict

#-----------------------------------------------------------------------
def read_viirs_m_band(**kwargs):
	'''
	Function to read VIIRS M-band file
	
	'''
	### import numpy as np
	### from netCDF4 import Dataset
	
	print(' - Calling read_viirs_m_band...\n')
	# path of the rad file
	radFile	= kwargs.get('radFile', None)
	# list of rad variable to be read
	radVars	= kwargs.get('radVars', None)
	# path of the geo file
	geoFile	= kwargs.get('geoFile', None)
	# list of geo variables to be read
	geoVars	= kwargs.get('geoVars', None)
	# logical switch to determine whether return the navigation data
	navigation	= kwargs.get('navigation', True)
	
	non_refl_list = ['M12', 'M13', 'M14', 'M15', 'M16']
	outDict = {}
	
	if radFile is not None:
		if len(radVars) > 0:
		
			print(' - Reading', radFile, '\n')
			
			ncid = Dataset(radFile, 'r')
			varList = list(ncid['/observation_data'].variables)

			# also read some global attributes of the rad-file.
			outDict['GRingPointLatitude'] = ncid.time_coverage_end
			outDict['GRingPointLongitude'] = ncid.time_coverage_start		
			
			ncid.set_auto_scale(False)			
			# netCDF will automatically apply the scale_factor and the add_offset
			# to the variable if its metadata contains these two.
			# However, it will also apply to the invalid value or scale the 
			# for some variable, the conversion equation is not the same
			# with the build in one...
# 			ncid.set_auto_mask(False)
		
			for key in radVars:
				if key in varList:
					if 'quality_flags' in key:
						outDict[key] = ncid["/observation_data/" + key][:]
					elif 'uncert_index' in key:
						ucidx = ncid["/observation_data/" + key][:].astype(float)
						scale_factor = ncid["/observation_data/" + key].scale_factor
						# valid_max    = ncid["/observation_data/" + key].valid_max
						# valid_min    = ncid["/observation_data/" + key].valid_min
						# use valid max and min to filter out invalid value
						# ucidx[np.where((ucidx > valid_max) | (ucidx < valid_min))] = np.nan		
						uc			 = 1 + scale_factor * ucidx**2
						uc			 = np.ma.filled(uc, fill_value = np.nan)	
						outDict[key] = uc
					else:
						if key in non_refl_list:
						
							DN           = ncid["/observation_data/" + key][:]			
							lut = ncid["/observation_data/" + key + '_brightness_temperature_lut'][:]
							scale_factor = ncid["/observation_data/" + key].scale_factor
							add_offset   = ncid["/observation_data/" + key].add_offset
							
							rad		 	 = DN.astype(float)  * scale_factor + add_offset
							rad		     = np.ma.filled(rad, fill_value = np.nan)
									
							outDict[key + '_rad'] = rad
							outDict['BT'+key ]    = lut[DN]
							
							outDict[key + '_LUT'] = lut[:]
							outDict['BT'+key ][np.where( outDict[key + '_rad'] != outDict[key + '_rad']) ] = np.nan						
							outDict['BT'+key ]     = np.ma.filled(outDict['BT'+key ], fill_value = np.nan)
					
						if key not in non_refl_list:
							refl = ncid["/observation_data/" + key][:].astype(float)
							scale_factor = ncid["/observation_data/" + key].scale_factor
							add_offset   = ncid["/observation_data/" + key].add_offset
							refl		 = refl * scale_factor + add_offset
							refl		 = np.ma.filled(refl, fill_value = np.nan)
							outDict[key + '_refl'] = refl							
						
							rad = ncid["/observation_data/" + key][:].astype(float)
							scale_factor = ncid["/observation_data/" + key].radiance_scale_factor
							add_offset   = ncid["/observation_data/" + key].radiance_add_offset
							rad		 	 = rad * scale_factor + add_offset
							rad		 	 = np.ma.filled(rad, fill_value = np.nan)
							outDict[key + '_rad'] = rad	
							
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')
			ncid.close()
	
	
	if geoFile is not None:
		if len(geoVars) > 0:
			ncid = Dataset(geoFile, 'r')
			geoList = list(ncid['/geolocation_data'].variables)
			print(' - Reading', geoFile, '\n')
			
			# also read some global attributes of the geo-file.
			outDict['GRingPointLatitude'] = ncid.GRingPointLatitude
			outDict['GRingPointLongitude'] = ncid.GRingPointLongitude
		
			outDict['SouthBoundingCoordinate'] = ncid.SouthBoundingCoordinate
			outDict['NorthBoundingCoordinate'] = ncid.NorthBoundingCoordinate
			outDict['WestBoundingCoordinate'] = ncid.WestBoundingCoordinate
			outDict['EastBoundingCoordinate'] = ncid.EastBoundingCoordinate

			outDict['DayNightFlag'] = ncid.DayNightFlag
			outDict['startDirection'] = ncid.startDirection
			outDict['endDirection'] = ncid.endDirection
			
			for key in geoVars:
				if key in geoList:
					temp = ncid["/geolocation_data/" + key][:]	
					outDict[key] = np.ma.filled(temp.astype(float), fill_value = np.nan)
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')

			if navigation:
				temp = ncid["/navigation_data/earth_moon_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_moon_distance'] = temp

				temp = ncid["/navigation_data/earth_sun_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_sun_distance'] = temp

			ncid.close()

	return outDict

#-----------------------------------------------------------------------	
def read_viirs_dnb(**kwargs):
	'''
	Function to read VIIRS M-band file
	
	'''
	### import numpy as np
	### from netCDF4 import Dataset
	
	print(' - Calling read_viirs_dnb...\n')
	# path of the rad file
	radFile	= kwargs.get('radFile', None)
	# list of rad variable to be read
	radVars	= kwargs.get('radVars', None)
	# path of the geo file
	geoFile	= kwargs.get('geoFile', None)
	# list of geo variables to be read
	geoVars	= kwargs.get('geoVars', None)
	# logical switch to determine whether return the navigation data
	navigation	= kwargs.get('navigation', True)
	
	outDict = {}
	
	if radFile is not None:
		if len(radVars) > 0:
			ncid = Dataset(radFile, 'r')
			obVarList = list(ncid['/observation_data'].variables)
			ncid.set_auto_mask(False)
			timeEnd = ncid.time_coverage_end
			timeStart = ncid.time_coverage_start	
			for key in radVars:
				if key in obVarList:
					if key =='DNB_observations':
						FillValue = ncid["/observation_data/" + key]._FillValue
						valid_max = ncid["/observation_data/" + key].valid_max
						valid_min = ncid["/observation_data/" + key].valid_min
						outDict[key] = ncid["/observation_data/" + key][:]			
						key_mask = key + '_mask'
						mask = np.ones_like(outDict[key])
						mask[np.where(outDict[key] == FillValue)] = 0
						mask[np.where((outDict[key] > valid_max) | (outDict[key] < valid_min))] = 0
						outDict[key_mask] = mask
					else:
						outDict[key] = ncid["/observation_data/" + key][:]	
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')
	
	
	if geoFile is not None:
		if len(geoVars) > 0:
			ncid = Dataset(geoFile, 'r')
			geoList = list(ncid['/geolocation_data'].variables)
			print(' - Reading', geoFile)
			
			# also read some global attributes of the geo-file.
			outDict['GRingPointLatitude'] = ncid.GRingPointLatitude
			outDict['GRingPointLongitude'] = ncid.GRingPointLongitude
		
			outDict['SouthBoundingCoordinate'] = ncid.SouthBoundingCoordinate
			outDict['NorthBoundingCoordinate'] = ncid.NorthBoundingCoordinate
			outDict['WestBoundingCoordinate'] = ncid.WestBoundingCoordinate
			outDict['EastBoundingCoordinate'] = ncid.EastBoundingCoordinate

			outDict['DayNightFlag'] = ncid.DayNightFlag
			outDict['startDirection'] = ncid.startDirection
			outDict['endDirection'] = ncid.endDirection
			
			outDict['OrbitNumber'] = ncid.OrbitNumber
			outDict['number_of_filled_scans'] = ncid.number_of_filled_scans
			
			outDict['PGE_StartTime'] = ncid.PGE_StartTime
			outDict['PGE_EndTime'] = ncid.PGE_EndTime
					
			outDict['platform'] = ncid.platform
			outDict['instrument'] = ncid.instrument			
			outDict['ShortName'] = ncid.ShortName
			
			print(' - DayNightFlag: '    , outDict['DayNightFlag'])
			print(' - Start Direction: ' , outDict['startDirection'])
			print(' - End   Direction: ' , outDict['endDirection'])
			print(' - Orbit Number: '    , outDict['OrbitNumber'])
			print(' - Start Time: '      , outDict['PGE_StartTime'])
			print(' - End Time: '        , outDict['PGE_EndTime'])
			print('\n')
			
			for key in geoVars:
				if key in geoList:
					temp = ncid["/geolocation_data/" + key][:]	
					outDict[key] = np.ma.filled(temp, fill_value = np.nan)
				else:
					print(' - Warning: Variable ' + key + ' is not in the dataset, skipping...')

			if navigation:
				temp = ncid["/navigation_data/earth_moon_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_moon_distance'] = temp

				temp = ncid["/navigation_data/earth_sun_distance"][:]
				temp = np.ma.filled(temp, fill_value = np.nan)
				temp = np.reshape(temp, (len(temp), 1))
				outDict['earth_sun_distance'] = temp	

			ncid.close()

	return outDict
#-----------------------------------------------------------------------
def read_VNP14(obFile, **kwargs):
	'''
	read nasa viirs VPN14, Thermal Anomalies/Fire 6-Min L2 Swath 750m
	'''
	### import numpy as np
	### from netCDF4 import Dataset	
	
	obVar	= kwargs.get('obVar', [])
	
	print(' - Reading ', obFile)
	
	output ={}
	
	hdfFile = Dataset(obFile)
	numFire = hdfFile.FirePix
	if numFire == 0:
		print(' - Warning: Not fire pixel in this file, return False...')
		return False
	
# 	print(hdfFile.variables)
	
	if len(obVar) > 0:
		obVarList = list(hdfFile.variables)
		hdfFile.set_auto_mask(False)
		for key in obVar:
			if key in obVarList:					
				output[key] = hdfFile[key][:]
			else:
				print(' - Warning read_VNP14: Variable '\
				      + key + ' is not in the dataset, skipping...')
		
	output['latitude'] = hdfFile['FP_latitude'][:]
	output['longitude'] = hdfFile['FP_longitude'][:]
	hdfFile.close()

	return output

#-----------------------------------------------------------------------
def read_Geos_FP(filename, GeosVar):
	'''
	read nasa Geos-FP surface flux data field
	'''

	### from netCDF4 import Dataset
	
	output = {}

	
	ncid = Dataset(filename)
	ncid.set_auto_mask(True)
	resol_lat = float(ncid.LatitudeResolution)
	resol_lon = float(ncid.LongitudeResolution)
	print(f' - Reading {filename}, resolution [Lat X Lon]: {resol_lat} X {resol_lon} \n')

	
	GeosList = list(ncid.variables)
	output['lat'] = ncid['lat'][:]
	output['lon'] = ncid['lon'][:]

	
	for key in GeosVar:
		if key in GeosList:
			output[key] = ncid[key][:]
		else:
			print('  Warning: Variable ' + key + ' is not in the dataset, skipping...')
			
	output['resol_lat'] = resol_lat
	output['resol_lon'] = resol_lon
	
	return output


#   :LatitudeResolution = "0.25";
#   :LongitudeResolution = "0.3125";

#   :LatitudeResolution = "0.5";
#   :LongitudeResolution = "0.625";

#-----------------------------------------------------------------------
def read_DayNightFlag(filename):

	### from netCDF4 import Dataset
	
	ncid = Dataset(filename, 'r')
	DayNightFlag = ncid.DayNightFlag

	ncid.close()
	
	
	return DayNightFlag

#-----------------------------------------------------------------------
def read_MCD12Q1(filename, params):
	
	'''
	Read the MCD12Q1 with pyhdf API
	
	
	'''
	
	from pyhdf.SD import SD, SDC
	import numpy as np
	
	print(' - read_MCD12Q1:', filename)
	
	# open the HDF4 file
	hdf = SD(filename, SDC.READ)
	
	output = {}
	
	for param in params:
		
		# get the Land_Cover_Type_1 data
		sds = hdf.select('LC_Type1')
		data = sds.get()
		
		# get the valid range and fill value metadata
		valid_range = sds.attributes()['valid_range']
		fill_value = sds.attributes()['_FillValue']

		# apply the valid range mask and fill value mask
		data = np.ma.masked_outside(data, valid_range[0], valid_range[1])
		data = np.ma.masked_equal(data, fill_value)	

		output[param] = data

	
	return output
#-----------------------------------------------------------------------
def read_nc(obFile, obVars):

	### from netCDF4 import Dataset
	print(' - Reading', obFile)
	ncid = Dataset(obFile)
	
	varList = list(ncid.variables)
	output  = {}

	for key in obVars:
		if key in varList:
			output[key] = ncid[key][:]
	ncid.close()
	return output	


#-----------------------------------------------------------------------	
def write_nc_mod(data_dict, savename):
	
	
	'''
	
	Retired, replaced by a more general approach write_nc and aux_infor_dict
	
	'''
	nFire = len(data_dict['FP_latitude'])
	
	save_dict = {}
	save_dict['FP_latitude']                = np.array( data_dict['FP_latitude'] )
	save_dict['FP_latitude_long_name']		= 'Latitude'
	save_dict['FP_latitude_unit'] 			= 'degree'

	save_dict['FP_longitude']               = np.array( data_dict['FP_longitude'] )
	save_dict['FP_longitude_long_name']     = 'Longitude'
	save_dict['FP_longitude_unit']			= 'degree'

	save_dict['FP_lines']                	= np.array( data_dict['FP_lines'] )
	save_dict['FP_lines_long_name']			= 'granule line of fire pixel'			
	save_dict['FP_lines_unit']				= 'None'

	save_dict['FP_samples']                	= np.array(data_dict['FP_samples'])
	save_dict['FP_samples_long_name']		= 'granule sample of fire pixel'
	save_dict['FP_samples_unit']			= 'None'

	save_dict['FP_T13(K)']                	= np.array(data_dict['FP_T13(K)'])
	save_dict['FP_T13(K)_long_name']		= 'Brightness Temperature at M13'
	save_dict['FP_T13(K)_unit']				= 'K'

	save_dict['FP_T15(K)']                	= np.array(data_dict['FP_T15(K)'])
	save_dict['FP_T15(K)_long_name']		= 'Brightness Temperature at M15'
	save_dict['FP_T15(K)_unit']				= 'K'

	
	save_dict['FP_Area(m2)']				=  np.array(data_dict['FP_Area(m2)'])
	save_dict['FP_Area(m2)_long_name']		= 'Pixel Area'
	save_dict['FP_Area(m2)_unit']			= 'm2'

	save_dict['FP_Visible_Energy_Fraction'] 		    = np.array(data_dict['FP_Visible_Energy_Fraction'])
	save_dict['FP_Visible_Energy_Fraction_long_name']	= 'Visible Energy Fraction'
	save_dict['FP_Visible_Energy_Fraction_unit']		= 'None'
	
	save_dict['FP_Visible_Energy(MW)'] 		    		= np.array(data_dict['FP_Visible_Energy(MW)'])
	save_dict['FP_Visible_Energy(MW)_long_name']		= 'Visible Energy'
	save_dict['FP_Visible_Energy(MW)_unit']				= 'MW'

	save_dict['FP_power(MW)']				= np.array(data_dict['FP_power(MW)'])
	save_dict['FP_power(MW)_long_name']		= 'fire radiative power'
	save_dict['FP_power(MW)_unit']			= 'MW'

	save_dict['FP_T13_std(K)']				= np.array(data_dict['FP_T13_std(K)'])
	save_dict['FP_T13_std(K)_long_name']	= 'background brightness temperature difference mean absolute deviation at M13'
	save_dict['FP_T13_std(K)_unit']			= 'K'
	
	save_dict['FP_Modified_Combustion_Efficiency']			 = np.array(data_dict['FP_Modified_Combustion_Efficiency'])
	save_dict['FP_Modified_Combustion_Efficiency_long_name'] = 'Modified Combustion Efficiency'
	save_dict['FP_Modified_Combustion_Efficiency_unit']		 = 'None'
	
	save_dict['FP_T04(K)_1']			 = np.array(data_dict['FP_T04(K)_1'])
	save_dict['FP_T04(K)_1_long_name']	= 'Brightness Temperature at I04, Pixel 1'
	save_dict['FP_T04(K)_1_unit']		= 'K'

	save_dict['FP_T04(K)_2']			 = np.array(data_dict['FP_T04(K)_2'])
	save_dict['FP_T04(K)_2_long_name'] 	= 'Brightness Temperature at I04, Pixel 2'
	save_dict['FP_T04(K)_2_unit']		= 'K'
	
	save_dict['FP_T04(K)_3']			 = np.array(data_dict['FP_T04(K)_3'])
	save_dict['FP_T04(K)_3_long_name']	= 'Brightness Temperature at I04, Pixel 3'
	save_dict['FP_T04(K)_3_unit']		= 'K'
	
	save_dict['FP_T04(K)_4']			 = np.array(data_dict['FP_T04(K)_4'])
	save_dict['FP_T04(K)_4_long_name'] 	= 'Brightness Temperature at I04, Pixel 4'
	save_dict['FP_T04(K)_4_unit']		= 'K'

	save_dict['FP_T05(K)_1']			 = np.array(data_dict['FP_T05(K)_1'])
	save_dict['FP_T05(K)_1_long_name']  = 'Brightness Temperature at I05, Pixel 1'
	save_dict['FP_T05(K)_1_unit']		= 'K'
	
	save_dict['FP_T05(K)_2']			 = np.array(data_dict['FP_T05(K)_2'])
	save_dict['FP_T05(K)_2_long_name']  = 'Brightness Temperature at I05, Pixel 2'
	save_dict['FP_T05(K)_2_unit']		= 'K'
	
	save_dict['FP_T05(K)_3']			 = np.array(data_dict['FP_T05(K)_3'])
	save_dict['FP_T05(K)_3_long_name']  = 'Brightness Temperature at I05, Pixel 3'
	save_dict['FP_T05(K)_3_unit']		= 'K'
	
	save_dict['FP_T05(K)_4']			 = np.array(data_dict['FP_T05(K)_4'])
	save_dict['FP_T05(K)_4_long_name']  = 'Brightness Temperature at I05, Pixel 4'
	save_dict['FP_T05(K)_4_unit']		= 'K'

	save_dict['Bowtie_Flag']			 = np.array(data_dict['bowtie_flag'])
	save_dict['Bowtie_Flag_long_name']  = 'Bowtie detection flag, 1 = detect bowtie, 0 = normal pixel'
	save_dict['Bowtie_Flag_unit']		= 'None'

	save_dict['FP_Land_Type']			= np.array(data_dict['FP_Land_Type'])
	save_dict['FP_Land_Type_long_name']	= 'MODIS land surface type'
	save_dict['FP_Land_Type_unit']		= 'None'

	ncid = Dataset(savename, 'w', format='NETCDF4' )
	ncid.createDimension('nFire', nFire)


	for key in save_dict.keys():
		if 'FP_Radiance_' in key:
			continue
	# 	print(type(save_dict[key]))
		dataType = type(save_dict[key])

		if dataType == np.ndarray:
			print(key, dataType)
			nDim = len(np.shape(save_dict[key]))
			print(nDim)
			if nDim == 1:
				print('Saving')
				tempInstance = ncid.createVariable( key,'f4', ('nFire') )
				tempInstance[:] = save_dict[key][:]
				if str(key) + '_units' in save_dict.keys():
					tempInstance.units = save_dict[str(key) + '_units']
		
				if str(key) + '_long_name' in save_dict.keys():
					tempInstance.long_name = save_dict[str(key) + '_long_name']		


	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	ncid.source = 'FILDA detecion on S-NPP M-Band'
	ncid.close()
	
	
	return

#-----------------------------------------------------------------------
def write_nc_img(data_dict, savename):

	nFire = len(data_dict['FP_latitude'])
	
	save_dict = {}
	save_dict['FP_latitude']                = np.array(data_dict['FP_latitude'])
	save_dict['FP_latitude_long_name']		= 'Longitude'
	save_dict['FP_latitude_unit'] 			= 'degree'

	save_dict['FP_longitude']               = np.array(data_dict['FP_longitude'])
	save_dict['FP_longitude_long_name']     = 'Longitude'
	save_dict['FP_longitude_unit']			= 'degree'

	save_dict['FP_lines']                	= np.array(data_dict['FP_lines'])
	save_dict['FP_lines_long_name']			= 'granule line of fire pixel'			
	save_dict['FP_lines_unit']				= 'None'

	save_dict['FP_samples']                	= np.array(data_dict['FP_samples'])
	save_dict['FP_samples_long_name']		= 'granule sample of fire pixel'
	save_dict['FP_samples_unit']			= 'None'

# 	save_dict['FP_T04(K)']                	= data_dict['FP_T04(K)']
# 	save_dict['FP_T04(K)_long_name']		= 'Brightness Temperature at I04'
# 	save_dict['FP_T04(K)_unit']				= 'K'
# 
# 	save_dict['FP_T05(K)']                	= data_dict['FP_T05(K)']
# 	save_dict['FP_T15(K)_long_name']		= 'Brightness Temperature at M15'
# 	save_dict['FP_T15(K)_unit']				= 'K'

	save_dict['FP_Area(m2)']				= np.array(data_dict['FP_Area(m2)'])
	save_dict['FP_Area(m2)_long_name']		= 'Pixel Area'
	save_dict['FP_Area(m2)_unit']			= 'm2'

	save_dict['FP_Visible_Energy_Fraction'] 		    = np.array(data_dict['FP_Visible_Energy_Fraction'])
	save_dict['FP_Visible_Energy_Fraction_long_name']	= 'Visible Energy Fraction'
	save_dict['FP_Visible_Energy_Fraction_unit']		= 'None'
	
	save_dict['FP_Visible_Energy(MW)'] 		    		= np.array(data_dict['FP_Visible_Energy(MW)'])
	save_dict['FP_Visible_Energy(MW)_long_name']		= 'Visible Energy'
	save_dict['FP_Visible_Energy(MW)_unit']				= 'MW'

	save_dict['FP_power(MW)']				= np.array(data_dict['FP_power(MW)'])
	save_dict['FP_power(MW)_long_name']		= 'fire radiative power'
	save_dict['FP_power(MW)_unit']			= 'MW'

# 	save_dict['FP_T13_std(K)']				= data_dict['FP_T13_std(K)']
# 	save_dict['FP_T13_std(K)_long_name']	= 'background brightness temperature difference mean absolute deviation at M13'
# 	save_dict['FP_T13_std(K)_unit']			= 'K'
	
	save_dict['FP_Modified_Combustion_Efficiency']			 = np.array(data_dict['FP_Modified_Combustion_Efficiency'])
	save_dict['FP_Modified_Combustion_Efficiency_long_name'] = 'Modified Combustion Efficiency'
	save_dict['FP_Modified_Combustion_Efficiency_unit']		 = 'None'
	
	save_dict['FP_T04(K)']			 	= np.array(data_dict['FP_T04(K)'])
	save_dict['FP_T04(K)_long_name']	= 'Brightness Temperature at I04'
	save_dict['FP_T04(K)_unit']			= 'K'


	save_dict['FP_T05(K)']			 = np.array(data_dict['FP_T05(K)'])
	save_dict['FP_T05(K)_long_name'] = 'Brightness Temperature at I05, Pixel 1'
	save_dict['FP_T05(K)_unit']		 = 'K'

	save_dict['Bowtie_Flag']			= np.array(data_dict['bowtie_flag'])
	save_dict['Bowtie_Flag_long_name']  = 'Bowtie detection flag, 1 = detect bowtie, 0 = normal pixel'
	save_dict['Bowtie_Flag_unit']		= 'None'

	save_dict['num_fire']			 	= np.array(data_dict['num_fire'])
	save_dict['num_fire_long_name']     = 'Number of I band fires within a same M band host'
	save_dict['num_fire_unit']		    = 'None'
	
	
	save_dict['FP_Land_Type']			= np.array(data_dict['FP_Land_Type'])
	save_dict['FP_Land_Type_long_name']	= 'MODIS land surface type'
	save_dict['FP_Land_Type_unit']		= 'None'
	
	ncid = Dataset(savename, 'w', format='NETCDF4' )
	ncid.createDimension('nFire', nFire)

	for key in save_dict.keys():
		if 'FP_Radiance_' in key:
			continue
	# 	print(type(save_dict[key]))
		dataType = type(save_dict[key])

		if dataType == np.ndarray:
			nDim = len(np.shape(save_dict[key]))
			if nDim == 1:
				print('Saving', key)
				tempInstance = ncid.createVariable( key,'f4', ('nFire') )
				tempInstance[:] = save_dict[key][:]
				if str(key) + '_units' in save_dict.keys():
					tempInstance.units = save_dict[str(key) + '_units']
		
				if str(key) + '_long_name' in save_dict.keys():
					tempInstance.long_name = save_dict[str(key) + '_long_name']		


	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	ncid.source = 'FILDA detecion on S-NPP M-Band'
	ncid.close()
	
	
	return	


#-----------------------------------------------------------------------
def aux_infor_dict():
	'''
	aux_infor_dict provides metadata for the fire detection output
	
	'''

	infor_dict = {}
	infor_dict['FP_Line'] 	= ['Granule line of fire pixel',  'None',  'None', 'i2']
	infor_dict['FP_Sample'] = ['Granule sample of fire pixel', 'None',  'None', 'i2']
	
	infor_dict['FP_Longitude'] = ['Longitude of fire pixel' , 'degree',  'None', 'f4']
	infor_dict['FP_Latitude']  = ['Latitude of fire pixel' , 'degree',  'None', 'f4']
	
	
	infor_dict['FP_M13_WinSize']	= ['M13 (4.05 um) spatial window size for FRP calculation' , 'None',  'None', 'i2']
	infor_dict['FP_M13_Rad_Mean'] 		= ['Background mean M13 (4.05 um) radiance for FRP calculation' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M13_Rad_MAD'] 		= ['Background median absolute deviation of M13 (4.05 um) radiance', 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M13_Rad'] 		= ['M13 (4.05 um) radiance for FRP calculation' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M13_Rad_Num']        = ['Number of background M13 (4.05 um) pixel for FRP calculation' , 'None',  'None', 'i2']
	infor_dict['FP_M13_BT']         = ['M13 (4.05 um) brightness temperature of fire pixel' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_M13_BT_Mean'] 	= ['M13 (4.05 um) brightness temperature of background' , 'Kelvin',  'None', 'f4']
	
	infor_dict['FP_M07_Rad'] 		= ['M07 (0.865 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M07_Rad_Mean']   = ['Background mean M07 (0.865 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M07_Rad_Num']        = ['Number of background M07 (0.865 um) pixel' , 'None',  'None', 'i2']
	
	infor_dict['FP_M08_Rad'] 		= ['M08 (1.24 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M08_Rad_Mean']   = ['Background mean M08 (1.24 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M08_Rad_Num']    = ['Number of background M08 (1.24 um) pixel' , 'None',  'None', 'i2']	

	infor_dict['FP_M10_Rad'] 		= ['M10 (1.61 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M10_Rad_Mean']   = ['Background mean M10 (1.61 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M10_Rad_Num']    = ['Number of background M10 (1.61 um) pixel' , 'None',  'None', 'i2']

	infor_dict['FP_M11_Rad'] 		= ['M11 (2.25 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M11_Rad_Mean']   = ['Background mean M11 (2.25 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M11_Rad_Num']    = ['Number of background M11 (2.25 um) pixel' , 'None',  'None', 'i2']
	
	# MZ Mar 20, 2024, add for GMAO output
	infor_dict['FP_M12_Rad'] 		= ['M12 (3.7 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M12_Rad_Mean']   = ['Background mean M11 (3.7 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M12_Rad_Num']    = ['Number of background M11 (3.7 um) pixel' , 'None',  'None', 'i2']
	
	infor_dict['FP_M14_Rad'] 		= ['M14 (8.55 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M14_Rad_Mean']   = ['Background mean M14 (8.55 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M14_Rad_Num']    = ['Number of background M14 (8.55 um) pixel' , 'None',  'None', 'i2']

	infor_dict['FP_M15_Rad'] 		= ['M15 (10.76 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M15_Rad_Mean']   = ['Background mean M15 (10.76 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M15_Rad_Num']    = ['Number of background M15 (10.76 um) pixel' , 'None',  'None', 'i2']

	infor_dict['FP_M16_Rad'] 		= ['M16 (12.00 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M16_Rad_Mean']   = ['Background mean M16 (12.00 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_M16_Rad_Num']    = ['Number of background M16 (12.00 um) pixel' , 'None',  'None', 'i2']
	
	# MZ Mar 20, 2024, add for GMAO output
	infor_dict['FP_I04_Rad']        = ['I04 (3.74 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_I04_Rad_Mean']   = ['Background mean I04 (3.74 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_I04_Rad_Num']    = ['Number of background I04 (3.74 um) pixel' , 'None',  'None', 'i2']
	# MZ Mar 20, 2024, add for GMAO output
	infor_dict['FP_I05_Rad']        = ['I05 (11.45 um) radiance of the fire pixel' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_I05_Rad_Mean']   = ['Background mean I05 (11.45 um) radiance' , 'W/m^2/sr/um',  'None', 'f4']
	infor_dict['FP_I05_Rad_Num']    = ['Number of background I05 (11.45 um) pixel' , 'None',  'None', 'i2']

	infor_dict['FP_I04_BT']  = ['I04 (3.74 um) brightness temperature of fire pixel' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_I05_BT']  = ['I05 (11.45 um) brightness temperature of fire pixel' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_IMG_BTD'] = ['Brightness temperature of fire pixel (I04 - I05)' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_DNB_POS'] = ['DNB abnormal possibility' , 'None',  'None', 'f4']
	
	infor_dict['FP_WinSize'] = ['background window size for fire detection' , 'None',  'None', 'i1']
	infor_dict['FP_Status']  = ['Fire status' , 'None',  'None', 'i1']	
	
	infor_dict['FP_I04_Mean'] = ['I04 (3.74 um) brightness temperature of background' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_I05_Mean'] = ['I05 (11.45 um) brightness temperature of background' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_BTD_Mean'] = ['Mean background brightness temperature difference (I04 - I05)' , 'Kelvin',  'None', 'f4']
	
	infor_dict['FP_I04_MAD']	= ['Background I04 (3.74 um) brightness temperature mean absolute deviation' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_I05_MAD']	= ['Background I05 (11.45 um) brightness temperature mean absolute deviation' , 'Kelvin',  'None', 'f4']	
	infor_dict['FP_BTD_MAD']	= ['Background brightness temperature difference mean absolute deviation (I04 - I05)' , 'Kelvin',  'None', 'f4']

	infor_dict['FP_CM'] = ['Cloud Mask' , 'None',  '0: Cloud, 1: Clean', 'i1']

	infor_dict['FP_Power'] = ['Fire radiative power' , 'MW',  'None', 'f4']
	infor_dict['FP_VEF'] = ['Visible energy fraction' , 'None',  'None', 'f4']
	infor_dict['FP_VE'] = ['Visible energy' , 'MW',  'None', 'f4']
	infor_dict['FP_MCE'] = ['Modified combustion efficiency' , 'None',  'None', 'f4']
	infor_dict['FP_Power_QA'] = ['Quality flag for fire radiative power calculation' , 'None',  'None', 'i1']
	infor_dict['FP_Area'] = ['Pixel area' , 'm^2',  'None', 'f4']
	
	infor_dict['FP_Num_Fire'] = ['Number of I band fires within a same M band host pixel' , 'None',  'None', 'i1']
	
	
	infor_dict['DNB_observations']  = ['DNB radiance of fire pixel' , 'nW/sr/m^2', 'None', 'f4']
	infor_dict['FP_Land_Type']   = ['MODIS Land cover product MCD12C1' , 'None',  'None', 'i1']
	infor_dict['FP_Peatland']    = ['Peatland flag, 0: Unknown, 1: Contain Peatland' , 'None',  'None', 'i1']
	infor_dict['FP_Peatfrac']    = ['Peatland fraction' , 'None', 'None', 'f4']
	infor_dict['FP_Gas_Flaring'] = ['Gas Flaring flag' , 'None',  '0: Unknown, 1: Gas flaring pixel', 'i1']

	infor_dict['FP_BG_Temp']     = ['Estimated background Temperature' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_Fire_Temp']   = ['Estimated fire Temperature', 'Kelvin',  'None', 'f4']
	infor_dict['FP_Fire_Frac']   = ['Estimated fire fraction', 'None',  'None', 'f4']
	infor_dict['FP_Opt_Status']  = ['Optimization status' , 'None',  'None', 'i1']


	infor_dict['FP_T04_1'] = ['Brightness temperature at I04, Pixel 1' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T04_2'] = ['Brightness temperature at I04, Pixel 2' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T04_3'] = ['Brightness temperature at I04, Pixel 3' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T04_4'] = ['Brightness temperature at I04, Pixel 4' , 'Kelvin',  'None', 'f4']
	
	infor_dict['FP_T05_1'] = ['Brightness temperature at I05, Pixel 1' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T05_2'] = ['Brightness temperature at I05, Pixel 2' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T05_3'] = ['Brightness temperature at I05, Pixel 3' , 'Kelvin',  'None', 'f4']
	infor_dict['FP_T05_4'] = ['Brightness temperature at I05, Pixel 4' , 'Kelvin',  'None', 'f4']
	
	infor_dict['FP_Bowtie'] = ['Fraction of bowtie effect' , 'None',  'None', 'f4']
	
	infor_dict['Sensor_Zenith']  = ['View zenith angle' , 'degree',  'None', 'f4']
	infor_dict['Sensor_Azimuth'] = ['View azimuth angle', 'degree',  'None', 'f4']
	
	infor_dict['Solar_Zenith']  = ['Solar zenith angle', 'degree',  'None', 'f4']
	
	infor_dict['Fire_mask'] = ['Fire_mask', 'None', '\n0 not-processed (non-zero QF)\n1 bowtie\n2 glint\n3 water\n4 clouds\n5 clear land\n6 unclassified fire pixel\n7 low confidence fire pixel\n8 nominal confidence fire pixel\n9 high confidence fire pixel\n', 'i1']
	
	infor_dict['Algorithm_QA'] = ['Algorithm QA', 'bit field', 'Refer to ATBD', 'u4']
	
	infor_dict['FP_AdjCloud'] = ['Number of adjacent cloud pixels', 'None', 'None', 'i1']
	infor_dict['FP_AdjWater'] = ['Number of adjacent water pixels', 'None', 'None', 'i1']
	
	infor_dict['FP_confidence'] = ['Detection confidence', 'None', 'None', 'i1']
	
	infor_dict['FP_SAA_flag'] = ['South Atlantic Anomaly, 0: Normal, 1: SAA', 'None', 'None', 'i1']
	

	return infor_dict	
	
	
#-----------------------------------------------------------------------
def write_nc(save_dict, sat, savename):

	'''
	write_nc writes the provided dictionary into *.nc files

	'''
	
	nFire = len(save_dict['FP_Latitude'])
	nRow = np.shape(save_dict['Fire_mask'])[0]
	nCol = np.shape(save_dict['Fire_mask'])[1]
	nRow_img, nCol_img = save_dict['Algorithm_QA'].shape
	
	if nFire <=0:
		print(' - FILDA: No fire detected...')
		return 
	
	
	infor_dict = aux_infor_dict()
	
	for key in save_dict.keys():
		save_dict[key] = np.array(save_dict[key])

			
	ncid = Dataset(savename, 'w', format='NETCDF4' )
	ncid.createDimension('nFire', nFire)
	ncid.createDimension('nRow', nRow)
	ncid.createDimension('nCol', nCol)	
	ncid.createDimension('nRow_QA', nRow_img)
	ncid.createDimension('nCol_QA', nCol_img)	
	
	if nCol == 3200:
		chunksize_x = 16
		chunksize_y = 3200
	else:
		chunksize_x = 32
		chunksize_y = 6400	
	
	for key in save_dict.keys():
		
		dataType = type(save_dict[key])
# 		print('type:', dataType, np.shape(save_dict[key]) )
		if dataType == np.ndarray:
			nDim = len(np.shape(save_dict[key]))
			if nDim == 1:
				if key in infor_dict.keys():
# 					print(key, save_dict[key][:])
					dtype = infor_dict[key][3]
					tempInstance = ncid.createVariable( key, dtype, ('nFire'), zlib=True, complevel = 4 )
					tempInstance[:] = save_dict[key][:]
					
					tempInstance.units = infor_dict[key][1]
					tempInstance.long_name = infor_dict[key][0]
					tempInstance.legend    = infor_dict[key][2]
					tempInstance.data_type = dtype	
# 				if key in infor_dict.keys():
# 					tempInstance.units = infor_dict[key][1]
# 					tempInstance.long_name = infor_dict[key][0]
			if nDim == 2:
				if key in infor_dict.keys():
					dtype = infor_dict[key][3]
					if key == 'Fire_mask':
						tempInstance = ncid.createVariable( key, dtype, ('nRow', 'nCol'), zlib=True, complevel = 9) #, shuffle=True, chunksizes = (chunksize_x, chunksize_y)
					if key == 'Algorithm_QA':
						tempInstance = ncid.createVariable( key, dtype, ('nRow_QA', 'nCol_QA'), zlib=True, complevel = 9) #, shuffle=True, chunksizes = (32, 6400)
					tempInstance[:, :] = save_dict[key][:, :]	
					
					tempInstance.units     = infor_dict[key][1]
					tempInstance.long_name = infor_dict[key][0]
					tempInstance.legend    = infor_dict[key][2]
					tempInstance.data_type = dtype
	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	ncid.source = 'FILDA detecion on ' + sat
	ncid.close()
	print(' - FILDA: Saving ' + savename)
	
	return	
#-----------------------------------------------------------------------
def decoding_algorithm_QA(QA):

	'''
	decoding_algorithm_QA translates the AF algorithm QA into meaningful
	information dictionary

	'''
	
	
	bit_infor    = {}
	bit_infor[0] = 'Channel I1 quality'
	bit_infor[1] = 'Channel I2 quality'
	bit_infor[2] = 'Channel I3 quality'
	bit_infor[3] = 'Channel I4 quality'
	bit_infor[4] = 'Channel I5 quality'
	bit_infor[5] = 'Geolocation data quality'
	bit_infor[6] = 'Channel M13 quality'
	bit_infor[7] = 'Unambiguous fire'
	bit_infor[8] = 'Background pixel'
	bit_infor[9] = 'Bright pixel rejection'
	bit_infor[10] = 'Candidate pixel'
	bit_infor[11] = 'Scene background '
	bit_infor[12] = 'Test 1'
	bit_infor[13] = 'Test 2'
	bit_infor[14] = 'Test 3'
	bit_infor[15] = 'Test 4'
	bit_infor[16] = 'Pixel saturation condition'
	bit_infor[17] = 'Glint condition'
	bit_infor[18] = 'Potential South Atlantic magnetic anomaly pixel'
	bit_infor[19] = 'Fire pixel over water'
	bit_infor[20] = 'Persistence test'
	bit_infor[21] = 'Persistence test'
	bit_infor[22] = 'Residual bow-tie pixel'
	
	QA_mask = {}
	
	
	for i in range(23):
		mask_temp = QA & int(2**i)
		mask_temp = mask_temp>>i
		QA_mask[i] = mask_temp
	
	
	return QA_mask, bit_infor
	
#-----------------------------------------------------------------------
def writeCSV(filename, dataDic, badvalue = 'N/A'):
    ### import pandas as pd
    '''
    writeCSV writes a dictionary into *.csv files

    '''
    
    
    print(' - Writing', filename)
    df = pd.DataFrame.from_dict(dataDic,'columns')
    df.to_csv(filename, index=False, na_rep = 'N/A')
    return df

#-----------------------------------------------------------------------   
def write_nc_NTL(savename, nData, **kwargs):
	import sys
	import time
	import numpy as np
	import netCDF4 as nc
	from netCDF4 import Dataset
	
	
	
	print(' - Writing', savename)
	source		= kwargs.get('source', '')
	string2attr	= kwargs.get('string2attr', True)
	complevel	= kwargs.get('complevel', 9)
	
	
	ncid = Dataset(savename, 'w', format='NETCDF4' )
	
	ncid.history = "Created " + time.ctime(time.time()) + ' ' + sys.argv[0]
	ncid.source = source

	nRow = np.shape(nData['latitude'])[0]
	nCol = np.shape(nData['latitude'])[1]
# 	nCol = np.shape(nData['longitude'])[0]
	
	ncid.createDimension('latitude', nRow)
	ncid.createDimension('longitude', nCol)

	tempInstance = ncid.createVariable( 'latitude','f4', ('latitude') )
	tempInstance[:] = nData['latitude'][:, 0]

	tempInstance = ncid.createVariable( 'longitude','f4', ('longitude') )
	tempInstance[:] = nData['longitude'][0, :]

	nData.pop('latitude')
	nData.pop('longitude')

	if 'level' in nData.keys():
		nLev = np.size(nData['level'])
		ncid.createDimension('lev', nLev)	

	for key in nData.keys():
		dataType = type( nData[key] )
		if dataType == np.ndarray:
			
			nDim = len(np.shape(nData[key]))
# 			print(np.shape(nData[key]))
			print(' - Dimention of dataset', key, ':' , nDim)
			if nDim == 3:
				tempInstance = ncid.createVariable( key,'f4', ('lev', 'lat', 'lon'), zlib=True,complevel = complevel)
				tempInstance[:, :, :] = nData[key][:, :, :]
				if str(key) + '_units' in nData.keys():
					tempInstance.units = nData[str(key) + '_units']
				
				if str(key) + '_long_name' in nData.keys():
					tempInstance.long_name = nData[str(key) + '_long_name']

			if nDim == 2:
				if np.shape(nData[key]) == (nRow, nCol):
					# usually we want to geo-locate the data...
					tempInstance = ncid.createVariable( key,'f4', ('latitude','longitude') ,zlib=True, complevel = complevel, chunksizes = (2400, 2400))
					tempInstance[:, :] = nData[key][:, :]
				
					if str(key) + '_units' in nData.keys():
						tempInstance.units = nData[str(key) + '_units']
					
					if str(key) + '_long_name' in nData.keys():
						tempInstance.long_name = nData[str(key) + '_long_name']
				else:
					# if the data can not be geo-located...
					# we thus define the data by itself...
					rowName = 'nRow_' + key
					colName = 'nCol_' + key
					ncid.createDimension(rowName, np.shape(nData[key])[0])
					ncid.createDimension(colName, np.shape(nData[key])[1])
					tempInstance = ncid.createVariable( key, 'i4', (rowName, colName), zlib=True,complevel = complevel)
					tempInstance[:, :] = nData[key][:, :]
				
					if str(key) + '_units' in nData.keys():
						tempInstance.units = nData[str(key) + '_units']
					
					if str(key) + '_long_name' in nData.keys():
						tempInstance.long_name = nData[str(key) + '_long_name']
			if nDim == 1:
				tempInstance = ncid.createVariable( key,'f4', ('lev'), zlib=True,complevel = complevel)
				tempInstance[:] = nData[key][:]
				if str(key) + '_units' in nData.keys():
					tempInstance.units = nData[str(key) + '_units']
				
				if str(key) + '_long_name' in nData.keys():
					tempInstance.long_name = nData[str(key) + '_long_name']			
			
			
		if dataType == str:
			if '_long_name' in key:
				continue
			if '_units' in key:
				continue
				
			if string2attr:
				# here we can  write string into  attributes...
				ncid.setncattr(key,  nData[key])
			else:
				# here we can aslo write string into data 
				ncid.createDimension(key, len(nData[key]))
				tempInstance = ncid.createVariable(key, 'S1', (key))
				tempInstance[:] = nc.stringtochar(np.array([nData[key]], 'S'))					


	ncid.close()
	return
	
#-----------------------------------------------------------------------
def copy_static(geo_data, namelist, time):

	
	print(f' - User request to copy static dataset to {namelist["RUN_DIR"]}')
	
	cord = [geo_data['NorthBoundingCoordinate'], geo_data['SouthBoundingCoordinate'], geo_data['WestBoundingCoordinate'], geo_data['EastBoundingCoordinate']]

	PlateCarree_tiles = FILDA_Time_Cord.get_tiles(cord)
	print(' - The following PlateCarree tiles will be copy to ', namelist['RUN_DIR'])
	print(' -', ' '.join(PlateCarree_tiles))

	Sinu_tiles, _, _, _, _ = FILDA_Time_Cord.get_tile_sinusoidal_2(cord, namelist, numCeil = 2400)
	print(' - The following Sinusoidal tiles will be copy to ', namelist['RUN_DIR'])
	print(' -', ' '.join(Sinu_tiles))
	
	
	#-------------------------------------------------------------------
	# get the filenames for the nighttime climatology
	NTL_DIR = namelist['NTL_DIR']
	NTL_PREFIX = namelist['NTL_prefix']
	if NTL_PREFIX == 'N/A':
		print(f' - get_peatland: No prefix for peatland is specified...')
		NTL_PREFIX = ''
			
# 	year  = int(float(time.Y))
# 	month = int(float(time.M))
	
	# make the year...
	if float(time.Y) >= namelist['Default_year_NTL']:
		NTL_year = str(int(namelist['Default_year_NTL']))
	else:
		NTL_year = time.Y

	NTL_JDN = FILDA_Time_Cord.JulianDay(NTL_year+time.M+'01', outtype = 'nasa')

	NTL_files = []
	for tile in PlateCarree_tiles:
# 		filename = glob.glob(NTL_DIR + '*' + NTL_year +  time.M + '.' + tile + '.nc')
		# MZ May-11-2023, modify to standard name for the NTL
		filename = glob.glob(NTL_DIR + '*' + NTL_PREFIX + '*' + NTL_JDN + '.' + tile + '*.nc')
		print(NTL_DIR + '*' + NTL_PREFIX + '*' + NTL_JDN + '*.*' + tile + '*.nc')
		print(filename)
		if len(filename)>0:
			NTL_files.append(filename[0])
		else:
			print(' - Warning! CANNOT FIND',  NTL_JDN + '.' + tile + '.nc for lighttime light')
	

	#------------------------------------
	# get the filenames for the Land surface type
	# the land surface type follow the laads layout
	SURF_DIR = namelist['SURF_DIR'] + time.Y + '/' + '001/'
	
	# make the year...
	if float(time.Y) >= namelist['Default_year_SURF']:
		SURF_year = str(int(namelist['Default_year_SURF']))
	else:
		SURF_year = time.Y
	
	SURF_files = []
	for tile in Sinu_tiles:
		filename = glob.glob(SURF_DIR + '*' + tile + '*.hdf')
		if len(filename)>0:
			SURF_files.append(filename[0])
		else:
			print(' - Warning! CANNOT FIND',  time.Y + time.M + '.' + tile + '.nc for land surface type')

	#--------
	# get the filenames for the peatland
	PEATLAND_DIR  = namelist['PEATLAND_DIR']
	PEATLAND_files = []
	for tile in PlateCarree_tiles:
		filename = glob.glob(PEATLAND_DIR + '*' + tile + '.nc')
		if len(filename)>0:
			PEATLAND_files.append(filename[0])
		else:
			print(' - Warning! CANNOT FIND' + tile + ' for Peatland')
			
	#--------
	# get the filenames for the gas flaring
	GASFLARE_DIR = namelist['GASFLARE_DIR']
	GASFLARE_files = [GASFLARE_DIR + 'gas_flaring_2012_2020.nc']
	
	#--------
	# get the filenames for the look up table
	LUT_DIR = LUT_DIR = namelist['LUT_DIR']	
	LUT_files = [ LUT_DIR + namelist['platform'] + '_DNB2IMG_Resampling_Lookup_Table.nc',
				  LUT_DIR + namelist['platform'] + '_DNB2MOD_Resampling_Lookup_Table.nc', 
				  LUT_DIR + namelist['platform'] + '_MOD2DNB_Resampling_Lookup_Table.nc',
				  LUT_DIR + namelist['platform'] + '_MOD_Pixelareas_Lookup_Table.nc', 
				  LUT_DIR + 'Infrared.csv']


	print(f' - FILDA_IO: Modify static path to {namelist["RUN_DIR"]}')
	
	sub_dirs = ['NTL_DIR', 'SURF_DIR', 'PEATLAND_DIR', 'GASFLARE_DIR', 'LUT_DIR']
	sub_files = [NTL_files, SURF_files, PEATLAND_files, GASFLARE_files, LUT_files]
	
	for sub_dir, sub_file in zip(sub_dirs, sub_files):
		namelist[sub_dir] = namelist['RUN_DIR'] + sub_dir.split('_')[0] + '/'
		if os.path.isdir(namelist[sub_dir]) == False:
			print(' - FILDA: Making directory', namelist[sub_dir])
			os.mkdir(namelist[sub_dir])
		
		exist_files = os.listdir(namelist[sub_dir])
		for filename in sub_file:
			temp_str = filename.split('/')[-1]
		
			if temp_str not in exist_files:
				shutil.copy2(filename, namelist[sub_dir])
				print(f' - FILDA_IO: Copy {filename} to {namelist[sub_dir]}')
			else:
				print(f' - FILDA_IO: {temp_str} already exits')
				
	return namelist
#-----------------------------------------------------------------------
def ini_output_dir(namelist, time):

	'''
	ini_output_dir creates the out put path for a given time

	'''
	#-- it is no necessary to send the manelist dictionary, just send the out_dir path

	img_save_dir = namelist['OUT_DIR'] + 'IMG/'
	if os.path.isdir(img_save_dir) == False:
		print(' - FILDA: Making directory', img_save_dir)
		os.mkdir(img_save_dir)

	img_save_dir = img_save_dir + time.Y + '/'
	if os.path.isdir(img_save_dir) == False:
		print(' - FILDA: Making directory', img_save_dir)
		os.mkdir(img_save_dir)

	img_save_dir = img_save_dir + time.DOY + '/'
	if os.path.isdir(img_save_dir) == False:
		print(' - FILDA: Making directory', img_save_dir)
		os.mkdir(img_save_dir)


	mod_save_dir = namelist['OUT_DIR'] + 'MOD/'
	if os.path.isdir(mod_save_dir) == False:
		print(' - FILDA: Making directory', mod_save_dir)
		os.mkdir(mod_save_dir)

	mod_save_dir = mod_save_dir + time.Y + '/'
	if os.path.isdir(mod_save_dir) == False:
		print(' - FILDA: Making directory', mod_save_dir)
		os.mkdir(mod_save_dir)

	mod_save_dir = mod_save_dir + time.DOY + '/'
	if os.path.isdir(mod_save_dir) == False:
		print(' - FILDA: Making directory', mod_save_dir)
		os.mkdir(mod_save_dir)


	return img_save_dir, mod_save_dir
