'''
This is the main library of Fire LIght Detection Algorithm Version 2 (FILDA-2)

Main developer  : Meng Zhou	   (MZ)
Co-developer    : Lorena Castro (LC)

Institution		: Atmospheric and Environmental Research Lab, the University of Iowa
				  

Dependence:
FILDA-2 dependents on the following python libraries:

#	libraries		version
1	os 				Default		Miscellaneous operating system interfaces
2	sys				Default		System-specific parameters and functions
3	copy			Default		Shallow and deep copy operations
4	datetime		Default		Basic date and time types
9	multiprocessing	Default		Process-based parallelism
5	numpy			1.19.1		https://numpy.org/
6	scipy			1.4.1		https://scipy.org/
7	pandas			1.1.3		https://pandas.pydata.org/		
8	netCDF4			1.5.1.2		https://unidata.github.io/netcdf4-python/


'''

#-----------------------------------------------------------------------
# inport the build libraries used in FILDA2
import os
import copy
from copy import deepcopy
import numpy as np
from netCDF4 import Dataset
import glob

# Import support module for the FILDA main...
import FILDA_IO			# module for reading and writing files
import FILDA_Time_Cord  # module for manipulating the time and coordinates
import FILDA_Resample 	# modele for physical resampling from DNB to M/I band
import FILDA_Cloud 		# modele for cloud masking
import FILDA_CLT	    # modele for generating the surface light climatology
import FILDA_BT

def init_detection(dirNameList, string_handler = ['DIR', 'prefix'], non_read_handler = '#'):
	'''
	init_detection reads namelist.input for FILDA
	
	-----------
	Parameters:
	num_str		: number of string variables. After this number, variables 
				  will be convert to float
	
	dirNameList	: directory for the namelist.input
	-----------
	Returns:
	namelist	: contain all the information in the namelist.input
	
	
	LC: 07-02-2022, add the a parameter dirNameList for better code structure
	
	'''
	import pandas as pd
	#import os
	from datetime import datetime
	print('************************************************************')
	print('')
	print('        Fire Light Detection Algorithm (FILDA) V1.0         ')
	print('')
	print('************************************************************\n')

	# get the system time
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print(" - FILDA: Current Time ", current_time)
	
	
	#filename = "namelist.input"
	filename = dirNameList + "namelist.input"
	namelist = {}
	with open(filename, 'r') as f:
		lines = f.readlines()
		for i, line in enumerate(lines):
			
			l = line.strip()
			print(' ', line.replace("\n", ""))
			if non_read_handler in l:
				continue
				
			name, value = l.split(':')
			name  = name.strip()
			value = value.strip()
			namelist[name]  = value
			
			combined_condition = True
			
			for handler in string_handler:

				if handler not in name:
					combined_condition = combined_condition & True
				else:
					combined_condition = combined_condition & False

			if combined_condition:
				namelist[name] = float(value)

	print('************************************************************\n')
	#------------------------------
	# make the necessary directory
	#------------------------------
	if os.path.isdir(namelist['OUT_DIR']) == False:
		os.mkdir(namelist['OUT_DIR'])
		
	# MZ Sept.30, 2022 remove for the operational code
# 	if os.path.isdir(namelist['FIG_DIR']) == False:
# 		os.mkdir(namelist['FIG_DIR'])
	
	return namelist

#-----------------------------------------------------------------------
def sel_candidates_mod(modData, namelist, time):
	
	'''
	sel_candidates_mod selects the fire candidates on the M-band resolution
	
	-----------
	Parameters:
	modData		: dictionary for holding the M-band data
	namelist	: namelist dictionary
	time		: time string to determine the time span for the 
				  surface light climatology
	-----------
	Returns:
	cdt_fire_mod : dictionary for holding the fire candidates detected on
				   M-band resolution
	DNB_infor	 : dictionary for holding the probability test results on 
				   DNB
	'''
	
	from scipy.stats import norm
	
	cdt_fire_mod = []
	
	#---------------------------
	# Generate DNB climatologies
	#---------------------------
	dnb_clt_var = ['DNB_At_Sensor_Radiance_500m_mean', 'DNB_At_Sensor_Radiance_500m_std']
	dnb_clt = FILDA_CLT.get_dnb_clt(modData, dnb_clt_var, time, namelist)
	
	# use 10 nW cm^-2 sr^-1 to filter out city...
	idx_city = np.where( dnb_clt['DNB_At_Sensor_Radiance_500m_mean'] > 10)
	
	# normalized the realtime observations based on the climatology
	dnb_stat = (modData['DNB_observations'] - dnb_clt['DNB_At_Sensor_Radiance_500m_mean'])/dnb_clt['DNB_At_Sensor_Radiance_500m_std']
	del dnb_clt
	
	# calculating the possibilities
	posDNB  = norm.cdf(dnb_stat)
	del dnb_stat
	
	# this step follows Polivka 2016, use the 1% lowest DNB radiance to define a local DNB threshold 
	thres_DNB = sorted( modData['DNB_observations'][modData['DNB_observations']==modData['DNB_observations']])
	valid_num = int(np.nansum(modData['DNB_observations']==modData['DNB_observations'])*0.01)
	thres_DNB = np.nanmean( thres_DNB[0:valid_num] )
	print(f' - FILDA2: 1% lowest DNB radiance {thres_DNB}')


	# put data into dictionary for other application
	DNB_infor = {}
	DNB_infor['posDNB']           = posDNB
	DNB_infor['idx_city']         = idx_city
	DNB_infor['DNB_observations'] = modData['DNB_observations']
	DNB_infor['thres_DNB_BG']     = thres_DNB
	DNB_infor['thres_DNB']        = thres_DNB

	if thres_DNB<4:
		thres_DNB = 4
		DNB_infor['thres_DNB'] = thres_DNB
		
	# DNB band selection, and also put them in to the candidate fires list #
	# MZ 03/10/2022, remove the criterion land_water_mask
	idx = np.where( (modData['longitude']==modData['longitude'])         & \
	                (posDNB>=namelist['thres_DNB'])                      & \
	                (modData['CM_2']>0)                                  & \
	                ((modData['BTD_MOD']>=2) | (modData['BTM13']>=290))  & \
	                (modData['DNB_observations']>thres_DNB)				 & \
	                (modData['land_water_mask'] ==1)  )
	                
	# get the fire candidates on M band resolution.
	FP_lines   = idx[0]
	FP_samples = idx[1]

	cdt_fire_mod = {}
	cdt_fire_mod['FP_line']   = np.concatenate( (FP_lines * 2,   FP_lines * 2 + 1, FP_lines * 2	,      FP_lines * 2 + 1) )
	cdt_fire_mod['FP_sample'] = np.concatenate( (FP_samples * 2, FP_samples * 2,   FP_samples * 2 + 1, FP_samples * 2 + 1))
	

	return cdt_fire_mod, DNB_infor

#-----------------------------------------------------------------------	
def sel_candidates_img(imgData, DNB_infor, namelist):
	'''
	sel_candidates_img selects the fire candidates on the I-band resolution
	
	-----------
	Parameters:
	imgData		: dictionary for holding the I-band data
	DNB_infor	: DNB possibility test results
	namelist	: namelist dictionary

	-----------
	Returns:
	cdt_fire_mod : dictionary for holding the fire candidates detected on
				   I-band resolution

	'''

	def get_spatial_static(data_array, dim):
		'''
		Internal function to calculate the mean and standard deviation 
		for a given array using FFT
		
		Parameter:
		data_array, input array for calculating the statistic
		dim, full dimension of the spatial window
		
		Return:
		mean, std, num
		
		'''
		kernel = np.ones((dim,dim))

		num_data = np.zeros_like(data_array)
		# set the valid value to 1
		num_data[data_array==data_array] = 1
		# set the nan value to zero
		data_array[data_array!=data_array] = 0
		
		# use FFT to do the quick convolve
		# total sum
		neighbor_sum     = signal.fftconvolve(data_array, kernel, mode='same')
		# number sum
		num_neighbor     = signal.fftconvolve(num_data, kernel, mode='same')
		# number sum
		neighbor_sum_sq  = signal.fftconvolve(data_array**2, kernel, mode='same')

		mean = neighbor_sum / num_neighbor

		std  = (neighbor_sum_sq/num_neighbor - mean**2)**0.5

		return mean, std, num_neighbor

	
	from scipy.signal import convolve2d
	from scipy import signal
	
	# Get the DNB probability test results for I Band
	thres_DNB              = DNB_infor['thres_DNB']	
	posDNB_img             = np.full_like(imgData['latitude'], np.nan)
	posDNB_img[0::2, 0::2] = DNB_infor['posDNB']
	posDNB_img[1::2, 0::2] = DNB_infor['posDNB']
	posDNB_img[0::2, 1::2] = DNB_infor['posDNB']
	posDNB_img[1::2, 1::2] = DNB_infor['posDNB']

	DNB_img             = np.full_like(imgData['latitude'], np.nan)
	DNB_img[0::2, 0::2] = DNB_infor['DNB_observations']
	DNB_img[1::2, 0::2] = DNB_infor['DNB_observations']
	DNB_img[0::2, 1::2] = DNB_infor['DNB_observations']
	DNB_img[1::2, 1::2] = DNB_infor['DNB_observations']
	
	
	# SET flag for absolute and background fire test
	abs_flag               = np.full_like(imgData['latitude'], 0).astype(int)
	bg_fire_flag           = np.full_like(imgData['latitude'], 0).astype(int)
	ir_fire_flag           = np.full_like(imgData['latitude'], 0).astype(int)	
	
	#-------------------------------------------------------------------
	# first pick out the absolute fire...
	# absolute fire has a temperature > 320
	# some saturate fire will not be identify by the normal criterion...
	abs_idx	= np.where( 
						( imgData['CM'] >0 ) & ( imgData['latitude'] == imgData['latitude'] )						  		  &\
						( ( (imgData['BTI04']>320.) & (imgData['I04_quality_flags']==0) ) 								      |\
						( (imgData['BTI04']>=355.)  & (imgData['I04_quality_flags']==4) & (imgData['I05_quality_flags']==0) ) |\
						( (imgData['BTD_IMG']<0.)   & (imgData['BTI05']>310.) & (imgData['I05_quality_flags']==0) )           |\
						( (imgData['BTI04']<=209.)  & (imgData['BTI05']>335.) )												  )\
					  )
	abs_flag[abs_idx] = 1
	print(' - FILDA: Find out', np.shape(abs_idx)[1], 'absolute fires')

	#-------------------------------------------------------------------
	# then use the rigid criterion to select potential background fire
	# MZ, 03/10/2022, remove the criterion land_water_mask
	bg_fire_idx = np.where(( imgData['BTI04']   >= 300)     & \
	                       ( imgData['BTD_IMG'] >= 10)  	 & \
	                       ( imgData['CM'] > 0)   			 & \
	                       ( imgData['latitude'] == imgData['latitude'] ) ) #& (imgData['land_water_mask']==1)
	bg_fire_flag[bg_fire_idx] = 1
	print(' - FILDA: Find out', np.shape(bg_fire_idx)[1], 'potential background fires')
	
	# Determined the dynamic threshold for selecting candidates
	BTI04   = copy.deepcopy(imgData['BTI04'])
	BTD_IMG = copy.deepcopy(imgData['BTD_IMG'])
	# First exclude the absolute fire pixel
	BTI04[abs_idx] = np.nan
	BTD_IMG[abs_idx] = np.nan
	# Then, exclude the previous identified cloud and water pixel
	BTI04[imgData['CM'] <1] = np.nan
	BTD_IMG[imgData['CM'] <1] = np.nan
	
	#-----------------------
	# MZ 07-04-2022, add the dynamical test for the BTD selection...
	dim = 501
	#-----------------------
	# form the ocean sampling	
	BTI04_OCEAN = copy.deepcopy(BTI04)
	BTI04_OCEAN[imgData['land_water_mask']==1] = np.nan	
	BTI04_OCEAN_BG, STD_OCEAN_BG, NUM_OCEAN_BG = get_spatial_static(BTI04_OCEAN, dim)

	#-----------------------
	# form the land sampling
	ocean_idx = imgData['land_water_mask']!=1
	BTI04_LAND = copy.deepcopy(BTI04)
	BTI04_LAND[ocean_idx] = np.nan
	BTI04_BG, STD_BG, NUM_BG    = get_spatial_static(BTI04_LAND, dim)

	BTI04_BG[ocean_idx] = BTI04_OCEAN_BG[ocean_idx]
	STD_BG[ocean_idx]   = STD_OCEAN_BG[ocean_idx]
	NUM_BG[ocean_idx]   = NUM_OCEAN_BG[ocean_idx]

	#-----------------------
	# form the ocean sampling	
	BTD_IMG_OCEAN = copy.deepcopy(BTD_IMG)
	BTD_IMG_OCEAN[imgData['land_water_mask']==1] = np.nan	
	BTD_IMG_OCEAN_BG, BTD_IMG_STD_OCEAN_BG, BTD_IMG_NUM_OCEAN_BG = get_spatial_static(BTD_IMG_OCEAN, dim)

	#-----------------------
	# form the land sampling
	ocean_idx = imgData['land_water_mask']!=1
	BTD_IMG_LAND = copy.deepcopy(BTD_IMG)
	BTD_IMG_LAND[ocean_idx] = np.nan
	BTD_IMG_BG, BTD_IMG_STD_BG, BTD_IMG_NUM_BG = get_spatial_static(BTD_IMG_LAND, dim)

	
	BTD_IMG_BG[ocean_idx]       = BTD_IMG_OCEAN_BG[ocean_idx]
	BTD_IMG_STD_BG[ocean_idx]   = BTD_IMG_STD_OCEAN_BG[ocean_idx]
	BTD_IMG_NUM_BG[ocean_idx]   = BTD_IMG_NUM_OCEAN_BG[ocean_idx]
	#------------------------
	# adding a constrain to the dynamical window
	idx           = np.where(NUM_BG<10)
	BTI04_BG[idx] = namelist['thres_BTI04']
	
	diff = BTI04 - BTI04_BG - STD_BG*2.5
	BT   = np.ones_like(BTI04)
	BT[np.where(diff<0)] = 0
	
	diff = BTD_IMG - BTD_IMG_BG - BTD_IMG_STD_BG*2.5
	BTD  = np.ones_like(BTI04)
	BTD[np.where(diff<0)] = 0

	# Find out the visible light anomaly
	# 03/10/2022, remove the criterion land_water_mask
	# cloud (bright and have large temperature difference) will mimic the fire behavior in the BTD and DNB over ocean surface, 
	# use a fix temperature to screen out it...
	# use (BTD>=1) to replace (imgData['BTD_IMG'] >= 7.5), MZ 07-04-2022
	vis_ano_idx = np.where(
	                        ( posDNB_img>namelist['thres_DNB'])            & \
	                        ( imgData['CM'] > 0) 			               & \
							( (BTD>=1) & ( BT>=1 ) & (DNB_img>thres_DNB) ) & \
							( imgData['latitude'] == imgData['latitude'] ) )
						    # ( imgData['land_water_mask']==1)    & \) 
	
	# then use the normal criterion to select fire candidates...
	# 03/10/2022, remove the criterion land_water_mask
	icandi_idx = np.where((( imgData['BTI04']   >= namelist['thres_BTI04'])     & \
	                       ( imgData['BTD_IMG'] >= namelist['thres_BTD_IMG']))  & \
	                       ( imgData['CM'] > 0)   								& \
	                       ( imgData['latitude'] == imgData['latitude'] ) ) #& (imgData['land_water_mask']==1)
	                       
	print(' - FILDA: Find out', np.shape(icandi_idx)[1], 'thermal light anomaly')


	FP_lines   = np.concatenate( (vis_ano_idx[0], icandi_idx[0], abs_idx[0] ) )
	FP_samples = np.concatenate( (vis_ano_idx[1], icandi_idx[1], abs_idx[1] ) )	
	fire_img = np.concatenate( ( np.expand_dims(FP_lines, axis = 1), np.expand_dims(FP_samples, axis = 1) ), axis = 1)
	fire_img = np.unique(fire_img, axis = 0)		
	print(' - FILDA: Concatenating fires into ' ,np.shape(fire_img)[0], ' I band candidates')
	
	ir_fire_flag[icandi_idx] = 1
	
	# Some exceptions for no fire is found...
	if np.shape(fire_img)[1] == 0:
	
		cdt_fire_img = {}
		cdt_fire_img['FP_line']          =  []
		cdt_fire_img['FP_sample']        =  []
		cdt_fire_img['FP_abs_img']       =  []
		cdt_fire_img['FP_bg_img']        =  []	
		cdt_fire_img['FP_Dnbg_I04']      =  []
		cdt_fire_img['FP_Dnbg_std_I04']  =  []		
	
		return cdt_fire_img
	
	FP_lines     = fire_img[:,0]
	FP_samples   = fire_img[:,1]
	
	cdt_fire_img = {}
	cdt_fire_img['FP_line']         = FP_lines
	cdt_fire_img['FP_sample']       = FP_samples
	cdt_fire_img['FP_abs_img']	    = abs_flag[ FP_lines,  FP_samples]
	cdt_fire_img['FP_bg_img']	    = bg_fire_flag[FP_lines,  FP_samples]
	cdt_fire_img['FP_ir_img']	    = ir_fire_flag[FP_lines,  FP_samples]
	cdt_fire_img['FP_Dnbg_I04']	    = BTI04_BG[FP_lines,  FP_samples]
	cdt_fire_img['FP_Dnbg_std_I04'] = STD_BG[FP_lines,  FP_samples]
		
	return cdt_fire_img

#-----------------------------------------------------------------------
def sel_candidates(modData, imgData, namelist, time):
	'''
	sel_candidates integrate results from sel_candidates_img and sel_candidates_mod
	
	-----------
	Parameters:
	modData		: dictionary for holding the M-band data
	imgData		: dictionary for holding the I-band data
	DNB_infor	: DNB possibility test results
	namelist	: namelist dictionary
	time		: time string to determine the time span for the 
				  surface light climatology
	-----------
	Returns:
	cdt_fire			: dictionary for holding the fire candidates
	DNB_infor['posDNB']	: dictionary for holding the probability test 
                          results on DNB
	'''
	
	# it is necessary to keep the detection on M band
	cdt_fire_mod, DNB_infor = sel_candidates_mod(modData, namelist, time)
	cdt_fire_img 		    = sel_candidates_img(imgData, DNB_infor, namelist)
	
	print(' - FILDA: Find', len(cdt_fire_mod['FP_line']), 'on DNB band...')
	print(' - FILDA: Find', len(cdt_fire_img['FP_line']), 'on I band...')
	
	# record the absolute fire index...	
	fire_abs   = np.full(np.shape(imgData['latitude']), 0).astype(int)
	abs_idx    = np.where(cdt_fire_img['FP_abs_img'] == 1)
	if np.size(abs_idx)>0:
		fire_abs[cdt_fire_img['FP_line'][abs_idx], cdt_fire_img['FP_sample'][abs_idx]] = 1

	# record the background high temperature index...	
	fire_bg   = np.full(np.shape(imgData['latitude']), 0).astype(int)
	bg_idx    = np.where(cdt_fire_img['FP_bg_img'] == 1)
	if np.size(bg_idx)>0:
		fire_bg[cdt_fire_img['FP_line'][bg_idx], cdt_fire_img['FP_sample'][bg_idx]] = 1

	# record the pure IR...	
	fire_ir   = np.full(np.shape(imgData['latitude']), 0).astype(int)
	ir_idx    = np.where(cdt_fire_img['FP_ir_img'] == 1)
	if np.size(ir_idx)>0:
		fire_ir[cdt_fire_img['FP_line'][ir_idx], cdt_fire_img['FP_sample'][ir_idx]] = 1
	
	# record the dynamical background threshold
	dynamical_bgi4   = np.full(np.shape(imgData['latitude']), np.nan)
	dynamical_bgi4[cdt_fire_img['FP_line'], cdt_fire_img['FP_sample']] = cdt_fire_img['FP_Dnbg_I04']
	
	# record the dynamical background threshold
	dynamical_bgi4_std   = np.full(np.shape(imgData['latitude']), np.nan)
	dynamical_bgi4_std[cdt_fire_img['FP_line'], cdt_fire_img['FP_sample']] = cdt_fire_img['FP_Dnbg_std_I04']

	# concatenate the Fire detections on I band & M band
	FP_line   = np.concatenate( (cdt_fire_mod['FP_line'],   cdt_fire_img['FP_line']) )
	FP_sample = np.concatenate( (cdt_fire_mod['FP_sample'], cdt_fire_img['FP_sample']) )
	
	cdt_fire_img = np.concatenate( ( np.expand_dims(FP_line, axis = 1), np.expand_dims(FP_sample, axis = 1) ), axis = 1)
	
	if len(cdt_fire_img) >0:
		cdt_fire_img =  np.unique(cdt_fire_img, axis = 0)

		print(' - FILDA: Find', np.shape(abs_idx)[1], 'absolute fires on I band...')
		print(' - FILDA: Find', np.shape(bg_idx)[1],  'potential high temperature fires on I band...')
		print(' - FILDA: Find', len(cdt_fire_img), 'in total...')

		cdt_fire_img   = np.array(cdt_fire_img).astype(int)
		FP_line_img    = cdt_fire_img[:,0]
		FP_sample_img  = cdt_fire_img[:,1]
	
		FP_line_mod    = cdt_fire_img[:,0]//2
		FP_sample_mod  = cdt_fire_img[:,1]//2	
		
		
# 		# MZ. add the twilight test here.
# 		DNB_observations = copy.deepcopy( modData['DNB_observations'] )
# 		# remove the fire light
# 		DNB_observations[FP_line_mod, FP_sample_mod] = np.nan
# 		# remove the city light
# 		DNB_observations[ DNB_infor['idx_city'] ] = np.nan
# 		# find the twilight zone
# 		idx_twilight = np.where( (modData['solar_zenith']>namelist['twilight_ang_min']) & \
# 		                         (modData['solar_zenith']<namelist['twilight_ang_max']) & \
# 		                         (DNB_observations == DNB_observations))
		
		# make a copy of the DNB radiance again
		DNB_observations = copy.deepcopy( modData['DNB_observations'] )
		# remove the impact of the background radiance	
		DNB_observations = DNB_observations - DNB_infor['thres_DNB_BG']
		
# 		# get the twilight DNB radiance
# 		print('----')
# 		print(np.shape(idx_twilight))
# 		if np.size(idx_twilight) > 10:
# 			rad_DNB_BG_tw = np.nanmean(DNB_observations[idx_twilight])
# 			print( f' - Find twilight pixel, average radiance is {rad_DNB_BG_tw}')
# 			print(DNB_observations[idx_twilight])
# 			DNB_observations[idx_twilight] = DNB_observations[idx_twilight] + DNB_infor['thres_DNB_BG'] - rad_DNB_BG_tw
# 			print(DNB_observations[idx_twilight])


		cdt_fire = {}
		cdt_fire['FP_line_img'] 	= FP_line_img
		cdt_fire['FP_sample_img']	= FP_sample_img
	
		cdt_fire['FP_line_mod']		= FP_line_mod
		cdt_fire['FP_sample_mod']	= FP_sample_mod
		cdt_fire['BTD_IMG']				= imgData['BTD_IMG'][FP_line_img, FP_sample_img]
		cdt_fire['BTI04']				= imgData['BTI04'][FP_line_img, FP_sample_img]
		cdt_fire['BTI05']				= imgData['BTI05'][FP_line_img, FP_sample_img]
		cdt_fire['CM_IMG']				= imgData['CM'][FP_line_img, FP_sample_img]
		cdt_fire['FP_latitude_img']		= imgData['latitude'][FP_line_img, FP_sample_img]
		cdt_fire['FP_longitude_img']	= imgData['longitude'][FP_line_img, FP_sample_img]
		cdt_fire['FP_land_water_mask']	= imgData['land_water_mask'][FP_line_img, FP_sample_img]
		cdt_fire['SAA_flag']			= imgData['SAA_flag'][FP_line_img, FP_sample_img]
		
		# some flag for QA
		cdt_fire['FP_abs_img']          = fire_abs[FP_line_img, FP_sample_img]
		cdt_fire['FP_bg_img']           = fire_bg[FP_line_img, FP_sample_img]
		cdt_fire['FP_ir_img']           = fire_ir[FP_line_img, FP_sample_img]
		
		cdt_fire['FP_Dnbg_I04']	        = dynamical_bgi4[FP_line_img, FP_sample_img]
		cdt_fire['FP_Dnbg_STD_I04']	    = dynamical_bgi4_std[FP_line_img, FP_sample_img]

		cdt_fire['CM_MOD']			 = modData['CM'][FP_line_mod, FP_sample_mod]
		cdt_fire['DNB_observations'] = DNB_observations[FP_line_mod, FP_sample_mod]
		cdt_fire['FP_latitude_mod']	 = modData['latitude'][FP_line_mod, FP_sample_mod]
		cdt_fire['FP_longitude_mod'] = modData['longitude'][FP_line_mod, FP_sample_mod]
		cdt_fire['BTM13']			 = modData['BTM13'][FP_line_mod, FP_sample_mod]
		cdt_fire['BTD_MOD']			 = modData['BTD_MOD'][FP_line_mod, FP_sample_mod]
	
		cdt_fire['FP_posDNB']		 = DNB_infor['posDNB'][FP_line_mod, FP_sample_mod]
		
	
	else:
		cdt_fire = False
	
	return cdt_fire, DNB_infor['posDNB']

#-----------------------------------------------------------------------
def get_BG_IMG(imgData, modData, cdt_fire, masking = True):

	'''
	get_BG_IMG generate the background temperature dataset by removing
	the cloudy/absolute pixels for land and ocean surface.
	
	For the land sufface, ocean pixel will be removed, and vice versa.

	-----------
	Parameters:
	imgData		: dictionary for holding the I-band data
	modData		: dictionary for holding the M-band data
	cdt_fire	: dictionary for holding the candidate fire

	-----------
	Returns:
	bg_BT	    : dictionary for holding background temperature fields

	'''

	# create the background, deepcopy is needed
	
# 	band = ['I04', 'I05', 'M13', 'M11']
	bg_BT = {}
	bg_BT['BTI04_LAND'] 	 = copy.deepcopy( imgData['BTI04'] )
	bg_BT['BTI05_LAND'] 	 = copy.deepcopy( imgData['BTI05'] )
	bg_BT['BTD_IMG_LAND'] 	 = copy.deepcopy( imgData['BTD_IMG'] )

	bg_BT['BTI04_OCEAN'] 	 = copy.deepcopy( imgData['BTI04'] )
	bg_BT['BTI05_OCEAN'] 	 = copy.deepcopy( imgData['BTI05'] )
	bg_BT['BTD_IMG_OCEAN'] 	 = copy.deepcopy( imgData['BTD_IMG'] )
	
	bg_BT['BTM13_LAND'] 	 = copy.deepcopy( modData['BTM13'] )
	bg_BT['BTM13_OCEAN'] 	 = copy.deepcopy( modData['BTM13'] )

	# MZ 29-03 2022, standard product only exclude cloud...
	#-------------------------------------------------------------------
	# for land dataset...
	# find cloudy and water pixel
	icandi_idx = np.where( ( imgData['CM'] < 1 )  | ( imgData['land_water_mask']!=1 ) )
	print(' - FILDA: Exclude', np.shape(icandi_idx)[1], 'cloudy pixel from land backgroud...') # cloudy pixels 
	# remove the watery and cloudy fire pixels
	bg_BT['BTI04_LAND'][icandi_idx] = np.nan
	bg_BT['BTI05_LAND'][icandi_idx] = np.nan
	bg_BT['BTD_IMG_LAND'][icandi_idx] = np.nan
	
	#-------------------------------------------------------------------
	# For ocean dataset....
	# find cloudy and land pixel
	icandi_idx = np.where( ( imgData['CM'] < 1 )  | ( imgData['land_water_mask']==1 ) )
	print(' - FILDA: Exclude', np.shape(icandi_idx)[1], 'cloudy pixels from ocean backgroud...') # cloudy pixels 
	# remove the watery and cloudy fire pixels
	bg_BT['BTI04_OCEAN'][icandi_idx] = np.nan
	bg_BT['BTI05_OCEAN'][icandi_idx] = np.nan
	bg_BT['BTD_IMG_OCEAN'][icandi_idx] = np.nan	

	#-------------------------------------------------------------------
	# find out the absolute and high potential background fire
	# MZ, there is a debating whether to remove all the possible fire or just the high temperature fire...
	# currently, we remove the absolute fire and high possible background fire pixels
	idx = np.where (    (cdt_fire['FP_abs_img'] == 1) |  (cdt_fire['FP_bg_img'] ==1) )
	# find the the lines and samples need to be removed
	# debating to what degree sample should be remove...
	line_remove   = cdt_fire['FP_line_img'][idx]
	sample_remove = cdt_fire['FP_sample_img'][idx]

	# remove the fire in the land dataset
	print(' - FILDA: Exclude', len(line_remove), ' potential fire pixels from backgroud...')
	# remove the all fire pixels in the land dataset
	bg_BT['BTI04_LAND'][line_remove, sample_remove] = np.nan
	bg_BT['BTI05_LAND'][line_remove, sample_remove] = np.nan
	bg_BT['BTD_IMG_LAND'][line_remove, sample_remove]= np.nan	
	# remove the fire in the ocean dataset
	bg_BT['BTI04_OCEAN'][line_remove, sample_remove] = np.nan
	bg_BT['BTI05_OCEAN'][line_remove, sample_remove] = np.nan
	bg_BT['BTD_IMG_OCEAN'][line_remove, sample_remove] = np.nan		
	
	#-------------------------------------------------------------------
	# also process the the M13 dataset...
	# for land dataset
	# find cloudy and water pixel
	mcandi_idx = np.where( ( modData['CM'] < 1 ) | ( modData['land_water_mask'] !=1 ) )
	bg_BT['BTM13_LAND'][mcandi_idx] = np.nan

	# for ocean dataset
	# find cloudy and land pixel
	mcandi_idx = np.where( ( modData['CM'] < 1 ) | ( modData['land_water_mask'] ==1 ) )
	bg_BT['BTM13_OCEAN'][mcandi_idx] = np.nan
	
	return bg_BT

#-----------------------------------------------------------------------
def get_BG_IMG_2(imgData, modData, cdt_fire, masking = True):
	'''
	Retired, but keep it for now...
	'''
	# create the background
	bg_BT = {}
	bg_BT['BTI04'] 	 = copy.deepcopy( imgData['BTI04'] )
	bg_BT['BTI05'] 	 = copy.deepcopy( imgData['BTI05'] )
	bg_BT['BTD_IMG'] = copy.deepcopy( imgData['BTD_IMG'] )

	if masking:
		# MZ 29-03 2022, standard product only exclude cloud to detect the gas flaring over ocean
		icandi_idx = np.where( ( imgData['CM'] < 1 ) ) # | ( imgData['land_water_mask']!=1 )
		print(' - FILDA: Exclude', np.shape(icandi_idx)[1], 'watery and from backgroud...') # cloudy pixels 
		# remove the watery and cloudy fire pixels
		bg_BT['BTI04'][icandi_idx] = np.nan
		bg_BT['BTI05'][icandi_idx] = np.nan
		bg_BT['BTD_IMG'][icandi_idx] = np.nan

		print(' - FILDA: Exclude', len(cdt_fire['FP_line_img']), 'all potential fire pixels from backgroud...')
		# remove the all potential fire pixels
		bg_BT['BTI04'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan
		bg_BT['BTI05'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan
		bg_BT['BTD_IMG'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan	
	
	else:

		bg_BT['BTI04'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan
		bg_BT['BTI05'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan
		bg_BT['BTD_IMG'][cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = np.nan			
		
	return bg_BT

#-----------------------------------------------------------------------
def get_BG_MOD(modData, fire_pixel, main_band = 'M13', ancillary_bands = ['M07', 'M08', 'M10', 'M11', 'M12', 'M14', 'M15', 'M16', 'I04', 'I05']):
	'''
	get_BG_MOD gets the M band background radiance for FRP calculating
	
	-----------
	Parameters:
	modData		: dictionary for holding the M-band data
	fire_pixel	: dictionary for holding the detected fires

	-----------
	Returns:
	rad_BG		: dictionary for holding background radiances on M-band(2D)
	fire_pixel  : dictionary for holding information of fire pixels
	'''
	
	fire_pix_mod = np.concatenate( ( np.expand_dims(fire_pixel['FP_line_mod'], axis = 1), 
	                                 np.expand_dims(fire_pixel['FP_sample_mod'], axis = 1) ), axis = 1)
	fire_pix_mod =  np.unique(fire_pix_mod, axis = 0)
	
		
	bands = [main_band] + ancillary_bands

	for band in bands:
		key = f'FP_{band}_Rad'
		fire_pixel[key] = modData[f'{band}_rad'][fire_pixel['FP_line_mod'], fire_pixel['FP_sample_mod']]
	
	# deep copy the rad13 and remove the fire and cloud pixel...
	# MZ 01-30-2022, other IR bands are also included...
	rad_BG = {}

	rad_BG['rad13_LAND']  = deepcopy(modData['M13_rad'])
	rad_BG['rad13_OCEAN'] = deepcopy(modData['M13_rad'])
	
	for band in ancillary_bands:
		rad_BG['rad' + band + '_LAND'] = deepcopy(modData[band + '_rad'])
		rad_BG['rad' + band + '_OCEAN'] = deepcopy(modData[band + '_rad'])
	
	# copy the DNB observation in to the background statistic
	rad_BG['rad' + 'DNB' + '_LAND'] = deepcopy(modData['DNB_observations'])
	rad_BG['rad' + 'DNB' + '_OCEAN'] = deepcopy(modData['DNB_observations'])
	
	# remove the fire pixel
	for band in rad_BG.keys():
		rad_BG[band][fire_pix_mod[:,0], fire_pix_mod[:,1]] = np.nan

	# remove the cloud pixel... 
	invalid_idx = np.where( (modData['CM'] < 1) )  # | (modData['land_water_mask']!=1)
	for band in rad_BG.keys():
		rad_BG[band][invalid_idx] = np.nan	

	# remove the land/water pixel
	land_idx = np.where( modData['land_water_mask']  == 1 )
	water_idx = np.where( modData['land_water_mask'] != 1 )
	
	# remove the land/ocean pixels
	for band in rad_BG.keys():
		if 'OCEAN' in band:
			rad_BG[band][land_idx] = np.nan	
		else:
			rad_BG[band][water_idx] = np.nan	
	
	return rad_BG, fire_pixel

#-----------------------------------------------------------------------
def cal_bgstat(bg_BT, FP_line, FP_sample, namelist, data_type = '_LAND', win_factor = 1, verbose = False):	
	'''
	cal_bgstat calculates the background statistics for one fire candidate

	-----------
	Parameters:
	bg_BT		: dictionary for holding the background brightness temperature fields
	FP_line		: lines number of the current fire pixel in the dataset
	FP_sample	: sample number of the current fire pixel in the dataset
	namelist	: namelist dictionary

	-----------
	Returns:
	The return of the function is conditional. If enough background pixels
	are found to form statistically significant results, it will return
	
		bg_stat: dictionary for holding statistical results
	
	otherwise, return
	
		False
	
	'''
	
	# get the threshold from the namelist
	thres_num_min	= namelist['thres_num']
	thres_frac		= namelist['thres_frac']
	half_win_ini	= namelist['half_win_ini']//2
	half_win_max	= namelist['half_win_max'] * win_factor
	win_step        = namelist['win_step']
	
	# get the max line and max sample of the dataset...
	max_line, max_sample = np.shape(bg_BT['BTI04_LAND'])
	
	#-------------------------------------------------------------------
	# translate the namelist values to the actual values used
	half_win    = half_win_ini
	thres_local = (2 * half_win + 1)**2 * thres_frac
	thres_num   = np.min((thres_num_min, thres_local))
	valid_num   = 0
	flag        = 0
	
	if verbose:
		print(' - half_win', half_win, 'thres_num', thres_num)
	
	# loop until there are enough pixels to form the statistics
	while (valid_num<thres_num):

		#--------------------------------------------
		FP_line_ini = FP_line - half_win
		FP_line_ini = int( np.max((FP_line_ini, 0)))
	
		FP_line_end = FP_line + half_win + 1
		FP_line_end = int( np.min((FP_line_end, max_line)))
	
		#--------------------------------------------
		FP_sample_ini = FP_sample - half_win
		FP_sample_ini = int( np.max((FP_sample_ini, 0)))
	
		FP_sample_end = FP_sample + half_win + 1
		FP_sample_end = int( np.min((FP_sample_end, max_sample)) )
	
		#--------------------------------------------
		local_I04 = bg_BT['BTI04' + data_type][FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
		
		#True is treated as 1 and False is treated as 0
		valid_num = np.nansum(local_I04==local_I04)
		
		# update the threshold for the next loop
		half_win = half_win + win_step
		thres_local = (2 * half_win + 1)**2 * thres_frac
		thres_num   = np.min((thres_num_min, thres_local))
		if verbose:
			print(' - In loop ', half_win, 'thres_num', thres_num, 'valid_num', valid_num)
		
		# the exception to break out the loop
		if half_win > half_win_max:
			flag = 1
			break
	
	# determine the status based on the loop flag
	if flag == 0:
		# have enough pixels to form the statistic
		
		local_I04 = bg_BT['BTI04' + data_type][FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
		local_I05 = bg_BT['BTI05' + data_type][FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
		local_BTD = bg_BT['BTD_IMG' + data_type][FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
	
		bg_stat = {}
		# first order moment
		bg_stat['FP_I04_Mean'] = np.nanmedian(local_I04)
		bg_stat['FP_I05_Mean'] = np.nanmedian(local_I05)
		bg_stat['FP_BTD_Mean'] = np.nanmedian(local_BTD)
		
		# second order moment
		bg_stat['FP_I04_MAD'] = np.nanmedian( abs(local_I04 - np.nanmedian(local_I04)) ) 
		bg_stat['FP_I05_MAD'] = np.nanmedian( abs(local_I05 - np.nanmedian(local_I05)) )         
		bg_stat['FP_BTD_MAD'] = np.nanmedian( abs(local_BTD - np.nanmedian(local_BTD)) )   
		
		bg_stat['FP_WinSize'] = (half_win-win_step)*2
		
		return bg_stat
	
	else:
		# do not have enough pixels to form the statistic
		return False

#-----------------------------------------------------------------------
def identify_one(cdt_fire, bg_BT, counter, namelist, verbose = False):
	'''
	identify_one determines the ultimate status for one fire candidate  

	-----------
	Parameters:
	cdt_fire		: dictionary for holding the candidates fires
	bg_BT			: dictionary for holding the background brightness temperature fields
	counter			: ID for the candidates fire
	namelist		: namelist dictionary

	-----------
	Returns:
	bg_stat			: dictionary for holding statistical results and detection results
	
	'''

	
	
	FP_line_img   = cdt_fire['FP_line_img'][counter]
	FP_sample_img = cdt_fire['FP_sample_img'][counter]
	FP_abs_img 	  = cdt_fire['FP_abs_img'][counter]
	
	if verbose:
		print( ' - FILDA:', FP_line_img, FP_sample_img, cdt_fire['FP_land_water_mask'][counter])
	
	I04           = cdt_fire['BTI04'][counter]
	I05           = cdt_fire['BTI05'][counter]
	BTD           = cdt_fire['BTD_IMG'][counter]
	posDNB        = cdt_fire['FP_posDNB'][counter]
	
	if cdt_fire['FP_land_water_mask'][counter] == 0:
		data_type = '_OCEAN'
		if verbose:
			print( ' - FILDA: using water dataset to generate background statistics...', cdt_fire['FP_land_water_mask'][counter])
	else:
		data_type = '_LAND'
		if verbose:
			print( ' - FILDA: using land dataset to generate background statistics...', cdt_fire['FP_land_water_mask'][counter])
			
	#-----------------------------------
	# absolute test...
	#-----------------------------------
	if FP_abs_img == 1:
		bg_stat = {}
		bg_stat['FP_Status']     = 1
		bg_stat['QA_FLAG']       = 1
		bg_stat['FP_I04_Mean']   = -999
		bg_stat['FP_I05_Mean']   = -999
		bg_stat['FP_BTD_Mean']   = -999
		bg_stat['FP_I04_MAD']    = -999
		bg_stat['FP_I05_MAD']    = -999
		bg_stat['FP_BTD_MAD']    = -999
		bg_stat['FP_WinSize']    = namelist['half_win_ini_FRP']
		return bg_stat
	
	
	QA_FLAG = 1	
	#-----------------------------------
	# Calculating background statistics
	#-----------------------------------
	bg_stat = cal_bgstat(bg_BT, FP_line_img, FP_sample_img, namelist, data_type=data_type, verbose = verbose)
	if bg_stat == False:
		if verbose:
			print(' - FILDA: WARNING, CANNOT FIND ENOUGH PIXEL IN THE CLEAN BACKGROUND DATA SET!')
			print(' - FILDA: Enlarge the window to select the background')

		if bg_stat == False:
		
			bg_stat = {}
			bg_stat['FP_Status']     = -999
			bg_stat['QA_FLAG']       = -999
			bg_stat['FP_I04_Mean']   = -999
			bg_stat['FP_I05_Mean']   = -999
			bg_stat['FP_BTD_Mean']   = -999
			bg_stat['FP_I04_MAD']    = -999
			bg_stat['FP_I05_MAD']    = -999
			bg_stat['FP_BTD_MAD']    = -999
			bg_stat['FP_WinSize']    = -999
	
			return bg_stat

	#-----------------------------------
	# Contexture test
	#-----------------------------------
	# determine the confidence level based on the DNB test
	# MZ 23-May-2024, modify criterion, only apply the relaxed threshold to land and non-SSA affected region
	if (posDNB>namelist['thres_DNB']) & (cdt_fire['FP_land_water_mask'][counter] != 0) & (cdt_fire['SAA_flag'][counter] != 1) :
		sig_lev   = 2.5
		delta_temp= 7.5
		fire_type = 2
	else:
		sig_lev   = 3
		delta_temp= 9
		fire_type = 3	
	
	# forming the contexture thresholds
	thres_BTD_Dyn	= bg_stat['FP_BTD_Mean'] + sig_lev*bg_stat['FP_BTD_MAD']
	thres_BTD_Abs	= bg_stat['FP_BTD_Mean'] + delta_temp
	thres_I04_Dyn   = bg_stat['FP_I04_Mean'] + sig_lev*bg_stat['FP_I04_MAD']
	
	bg_stat['QA_FLAG'] = QA_FLAG
	
	flag_test_1 = 0
	if (BTD >= thres_BTD_Dyn):
		flag_test_1 = 1

	flag_test_2 = 0
	if (BTD >= thres_BTD_Abs):
		flag_test_2 = 1	

	flag_test_3 = 0
	if (I04 >= thres_I04_Dyn):
		flag_test_3 = 1
	
	# 
	total_flag = flag_test_1*2**2 + flag_test_2*2 + flag_test_3

	if total_flag == 7:
		bg_stat['FP_Status'] = total_flag
	else:
		bg_stat['FP_Status'] = -total_flag
	
	return bg_stat

#-----------------------------------------------------------------------
def MP_fire_test(cdt_fire, bg_BT, namelist, interest_pos):

	'''
	MP_fire_test is the basic detection unit for the multiprocessing...

	-----------
	Parameters:
	cdt_fire		: dictionary for holding the candidates fires
	bg_BT			: dictionary for holding the background brightness temperature fields
	counter			: ID for the candidates fire
	interest_pos    : section of the array that needs to be process in this call

	-----------
	Returns:
	result			: dictionary for holding the test results
	
	'''
	
	
	rows = interest_pos[1] - interest_pos[0]
	result = dict(x = [], FP_Status = [], QA_FLAG = [], FP_I04_Mean = [], \
				  FP_I05_Mean = [], FP_BTD_Mean = [], FP_I04_MAD = [], FP_I05_MAD = [], \
				  FP_BTD_MAD =[], FP_WinSize = [])


	for fire_counter in range(interest_pos[0], interest_pos[1]):
	
		one_result = identify_one(cdt_fire, bg_BT, fire_counter, namelist)

		result['x'].append(fire_counter)
		
		for var in one_result.keys():
			result[var].append(one_result[var])

	return result
	
#-----------------------------------------------------------------------
def fire_test(cdt_fire, bg_BT, namelist):

	'''
	fire_test do the fire test in a multiprocessing manner

	-----------
	Parameters:
	cdt_fire		: dictionary for holding the candidates fires
	bg_BT			: dictionary for holding the background brightness temperature fields

	-----------
	Returns:
	fire_pixel		: dictionary for holding the detected fire  pixels
	cdt_fire		: dictionary for holding the candidates fire pixels
	
	'''
	
	# import the necessary library here
	import multiprocessing
	from multiprocessing import Pool
	from multiprocessing import Process
	from functools import partial
	
	def split_array(len_array, processes_n):
		'''
		Internal function to split the job into multiple sub-jobs for 
		multiprocessing	
		'''
		list_limits = []
		for ii in np.arange(0, processes_n):
			total_items = len(len_array)
			q = int(total_items/processes_n)
			r = total_items % processes_n
			limits = []
	
			for iter in np.arange(0, processes_n-1):
				if (iter <= r-1) :
					limits.append(q+1)
				else: 
					limits.append(q)
			limits.append(total_items-sum(limits))
			lio = 0
			ix=ii
	
			while ( ix >= 1 ):
				lio = lio + limits[ix-1]
				ix = ix-1
			ls = lio + limits[ii] 
			lis = len_array[lio:ls]  
			list_limits.append([lio, ls])
		return list_limits	


	#-------------------------------------------------------------------
	# main code begin here
	num_fire = len(cdt_fire['FP_line_img'])
	
	if num_fire == 0:
		return False, False
	
	len_array = np.arange(0,num_fire)
	numCore = multiprocessing.cpu_count()
	
	# one shall modify here based on their super computer ability
	if num_fire < 500:
		numCore = 1
	elif (num_fire >= 500) & (num_fire < 5000):
		numCore = 5
	elif (num_fire >= 5000) & (num_fire < 20000):
		numCore = 10
	else:
		numCore = 20
		
	# create the array for the detection results
	FP_Status   = np.full( (num_fire, ), np.nan)
	QA_FLAG 	= np.full( (num_fire, ), np.nan)
	I04_mean	= np.full( (num_fire, ), np.nan)
	I05_mean	= np.full( (num_fire, ), np.nan)
	BTD_mean	= np.full( (num_fire, ), np.nan)
	I04_MAD 	= np.full( (num_fire, ), np.nan)
	I05_MAD		= np.full( (num_fire, ), np.nan)
	BTD_MAD		= np.full( (num_fire, ), np.nan)
	FP_WinSize	= np.full( (num_fire, ), np.nan)

	#------------------------------------------------
	# use multiprocessing to accelerate the speed...
	#------------------------------------------------
	sub_arrays_list = split_array(len_array, numCore)
	print(' - FILDA: Divide the fire candidates into', len(sub_arrays_list), 'sigments...')
	# set up the worker pools
	pool = Pool(processes=numCore)
	# define the partial functions for parallel computing
	func = partial(MP_fire_test, cdt_fire, bg_BT, namelist)
	# pool all the results together...
	result = pool.map(func, sub_arrays_list)
	# then concatenate all the parallel outcome...
	for item in result:
		FP_Status[item['x']] 	= item['FP_Status']
		QA_FLAG[item['x']] 	   	= item['QA_FLAG']
		I04_mean[item['x']] 	= item['FP_I04_Mean']
		I05_mean[item['x']]		= item['FP_I05_Mean']
		BTD_mean[item['x']]		= item['FP_BTD_Mean']
		I04_MAD[item['x']]		= item['FP_I04_MAD']
		I05_MAD[item['x']]		= item['FP_I05_MAD']
		BTD_MAD[item['x']]		= item['FP_BTD_MAD']
		FP_WinSize[item['x']]	= item['FP_WinSize']

	# close the pooling process...
	pool.close()
	pool.join()	
	
	# add the detection results to all the candidate fires
	cdt_fire['FP_Status']   = FP_Status
	cdt_fire['QA_FLAG']     = QA_FLAG
	cdt_fire['FP_I04_Mean'] = I04_mean
	cdt_fire['FP_I05_Mean'] = I05_mean
	cdt_fire['FP_BTD_Mean'] = BTD_mean
	cdt_fire['FP_I04_MAD']  = I04_MAD
	cdt_fire['FP_I05_MAD']  = I05_MAD
	cdt_fire['FP_BTD_MAD']  = BTD_MAD
	cdt_fire['FP_WinSize']  = FP_WinSize
	
	# forming the detected fire dataset
	# MZ 23-May-2024
	# One more QA constraint put here, M band brightness temperature should at least
	# be positive
	valid_idx = np.where( (cdt_fire['FP_Status']>0) & (cdt_fire['BTD_MOD']>0))
	
	fire_pixel = {}
	for key in cdt_fire.keys():
		fire_pixel[key] = cdt_fire[key][valid_idx]
	
	
	# Check the southern Atlantic anomaly here
	# MZ, Jun 01 2024, move this check to get_fire_paras
	# fire_pixel['FP_SAA_flag'] = check_saa(fire_pixel, bg_BT)
	# print(' - FILDA: Identify ', len(fire_pixel['FP_line_mod']), 'fires...')

	return fire_pixel, cdt_fire

#-----------------------------------------------------------------------
def check_saa(fire_pixel, bg_BT):

	'''
	check_saa checks whether a detected fire is caused by the southern 
	Atlantic anomaly

	-----------
	Parameters:
	fire_pixel		: dictionary for holding the detected fire  pixels
	bg_BT			: dictionary for holding the background brightness temperature fields

	-----------
	Returns:
	FP_SAA_CHECK    : checking status for each fire pixel
	'''

	# geht the M-band index for each fire
	line_remove = fire_pixel['FP_line_mod'].astype(int)
	sample_remove = fire_pixel['FP_sample_mod'].astype(int)

	bg_BT['BTM13_LAND'][line_remove, sample_remove] = np.nan
	bg_BT['BTM13_OCEAN'][line_remove, sample_remove] = np.nan

	print( ' - FILDA, checking SAA')
	FP_SAA_CHECK = []
	
	land_water_suffix = ['_OCEAN', '_LAND']
	
	# get the fire general surface type for each fire
	FP_data_types = [ land_water_suffix[idx] for idx in fire_pixel['FP_land_water_mask'].astype(int)]
	
	# at last, do an SAA check...
	for i in range(len(fire_pixel['SAA_flag'])):

		if fire_pixel['SAA_flag'][i] == 1:

			FP_line   = fire_pixel['FP_line_mod'][i]
			FP_sample = fire_pixel['FP_sample_mod'][i]
			data_type = FP_data_types[i]
			half_win  = int(fire_pixel['FP_WinSize'][i]//2 + 1)
			max_line, max_sample = np.shape(bg_BT['BTM13_OCEAN'])
		
			#--------------------------------------------
			FP_line_ini = FP_line - half_win
			FP_line_ini = int( np.max((FP_line_ini, 0)))

			FP_line_end = FP_line + half_win + 1
			FP_line_end = int( np.min((FP_line_end, max_line)))

			#--------------------------------------------
			FP_sample_ini = FP_sample - half_win
			FP_sample_ini = int( np.max((FP_sample_ini, 0)))

			FP_sample_end = FP_sample + half_win + 1
			FP_sample_end = int( np.min((FP_sample_end, max_sample)) )	
		
			line   = fire_pixel['FP_line_mod'][i]
			sample = fire_pixel['FP_sample_mod'][i]

			# create the threshold for M13 check...
			thres = np.nanmax(bg_BT['BTM13' + data_type][ FP_line_ini:FP_line_end,  FP_sample_ini:FP_sample_end]) + 2

			if thres == thres:		
				if fire_pixel['BTM13'][i] > thres:
					FP_SAA_CHECK.append(0)
				else:

					FP_SAA_CHECK.append(1)
			else:
				FP_SAA_CHECK.append(1)

		else:
		
			FP_SAA_CHECK.append(0)
	
	FP_SAA_CHECK = np.array(FP_SAA_CHECK)
	print(' - FILDA, number of Normal fire', np.size(FP_SAA_CHECK[np.where(FP_SAA_CHECK == 0)]))
	print(' - FILDA, number of SAA affected fire',np.size(FP_SAA_CHECK[np.where(FP_SAA_CHECK == 1)]))
	
	return FP_SAA_CHECK
	
	
#-----------------------------------------------------------------------
def get_fire_rad13(fire_pixel, rad_BG, namelist, modData):
	
	'''
	get_fire_rad13 gets the M13 (4.0 µm) get the background radiance for M13 to 
	calculate FRP on M band for I band detection

	-----------
	Parameters:
	rad_BG		: dictionary for holding background radiances on M-band(2D)
	fire_pixel  : dictionary for holding information of fire pixels
	namelist    : namelist dictionary
	modData		: dictionary for holding the M-band data
	-----------
	Returns:
	fire_pixel	: dictionary for holding information of the detected fire pixel
	'''

	def get_local_rad13(rad_BG, FP_line, FP_sample, FP_WinSize, data_type, namelist, modData):
		
		'''
		Internal function to get the local M13 (4.0 µm) radiance

		-----------
		Parameters:
		rad_BG		: dictionary for holding background radiances on M-band(2D)
		FP_line		: lines number of the current fire pixel in the dataset
		FP_sample	: sample number of the current fire pixel in the dataset
		FP_WinSize	: window size used to acquire the sample pixels
		data_type	: land/ocean string.
		namelist	: namelist dictionary
		modData		: dictionary for holding the M-band data
		-----------
		Returns:
		bg_rad_M13	: dictionary for holding the background  M13 (4.0 µm) radiance 
					  for each detected fire

		'''
		
		rad13           = rad_BG['rad13' + data_type]
		thres_num_min	= namelist['thres_num_FRP']
		thres_frac		= namelist['thres_frac_FRP']
		half_win_ini	= int(FP_WinSize//2 + 1)
		half_win_max	= namelist['half_win_max_FRP']
		win_step        = namelist['win_step_FRP']
		
		# get the max line and max sample of the dataset...
		max_line, max_sample = np.shape(rad13)
		#--------------------------------------------
		half_win    = half_win_ini
		thres_local = (2 * half_win + 1)**2 * thres_frac
		thres_num   = np.min((thres_num_min, thres_local))
		valid_num   = 0
		flag        = 0
		ancillary_record = {}
		ancillary_bands = ['M07', 'M08', 'M10', 'M11','M12', 'M14', 'M15', 'M16','I04', 'I05', 'DNB']

		while (valid_num<thres_num):
		
			#--------------------------------------------
			FP_line_ini = FP_line - half_win
			FP_line_ini = int( np.max((FP_line_ini, 0)))
	
			FP_line_end = FP_line + half_win + 1
			FP_line_end = int( np.min((FP_line_end, max_line)))
			
			#--------------------------------------------
			FP_sample_ini = FP_sample - half_win
			FP_sample_ini = int( np.max((FP_sample_ini, 0)))
	
			FP_sample_end = FP_sample + half_win + 1
			FP_sample_end = int( np.min((FP_sample_end, max_sample)) )
	
			#--------------------------------------------
			local_M13 = rad13[FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
			
			# True is treated as 1 and False is treated as 0 
			# use the lowest 25% of pixel as the the background
			valid_num = int(np.nansum(local_M13==local_M13)*0.25)
		
			# update the threshold for the next loop
			half_win = half_win + win_step
			thres_local = (2 * half_win + 1)**2 * thres_frac
			thres_num   = np.min((thres_num_min, thres_local))

		
			if half_win > half_win_max:
				flag = 1
				bg_rad_M13 = {}
				bg_rad_M13['FP_M13_Rad_Mean'] 	 = np.nan
				bg_rad_M13['FP_M13_Rad_MAD'] 	 = np.nan
				bg_rad_M13['FP_M13_Rad_Num'] 	 = np.nan
				bg_rad_M13['FP_M13_WinSize'] = half_win - 1
				bg_rad_M13['FP_Power_QA'] 	 = -1
				
				for ab in ancillary_bands:			
					bg_rad_M13['FP_' + ab + '_Rad_Mean'] = np.nan
					bg_rad_M13['FP_' + ab + '_Rad_Num']  = np.nan
					bg_rad_M13['FP_' + ab + '_Rad_STD']  = np.nan

				return bg_rad_M13

		
		#---------------------------------------------------------------
		# main code start here
		local_M13 = sorted(local_M13[local_M13==local_M13])[0:valid_num]
		bg_M13_mean = np.nanmean(local_M13)
		bg_M13_MAD  = np.nanmedian( abs(local_M13 - np.nanmedian(local_M13)) ) 


		for ab in ancillary_bands:
		
			ancillary_record[ab] = []
			
			local_ab =rad_BG['rad' + ab + data_type]
	
			local_ab = local_ab[FP_line_ini : FP_line_end, FP_sample_ini:FP_sample_end]
			valid_ab   = np.size(local_ab)

			if valid_ab > valid_num:
				valid_ab = valid_num
				local_ab   = sorted(local_ab[local_ab==local_ab])[0:valid_ab]
				bg_ab_mean = np.nanmean(local_ab)
				bg_ab_std  = np.nanstd(local_ab)
			elif (valid_ab > 0) & (valid_ab < valid_num):
				local_ab   = sorted(local_ab[local_ab==local_ab])[0:valid_ab]
				bg_ab_mean = np.nanmean(local_ab)
				bg_ab_std = np.nanstd(local_ab)
			else:
				valid_ab = np.nan
				bg_ab_mean = np.nan
				bg_ab_std = np.nan

			ancillary_record[ab] = [bg_ab_mean, bg_ab_std, valid_ab]

		half_win  = half_win - 1
		QA        = 1
		
		if half_win>31:
			QA = 0

		bg_rad_M13 = {}
		bg_rad_M13['FP_M13_Rad_Mean'] 	 = bg_M13_mean
		bg_rad_M13['FP_M13_Rad_MAD'] 	 = bg_M13_MAD
		bg_rad_M13['FP_M13_Rad_Num'] 	 = valid_num
		bg_rad_M13['FP_M13_WinSize'] = half_win
		bg_rad_M13['FP_Power_QA'] 	 = QA
		
		
		for key in ancillary_record.keys():
			
			bg_rad_M13['FP_' + key + '_Rad_Mean'] = ancillary_record[key][0]
			bg_rad_M13['FP_' + key + '_Rad_STD']  = ancillary_record[key][1]
			bg_rad_M13['FP_' + key + '_Rad_Num']  = ancillary_record[key][2]
		
		
		return  bg_rad_M13
	
	#---------------------------------------------------
	# main code start here
	FP_lines	= fire_pixel['FP_line_mod']
	FP_samples	= fire_pixel['FP_sample_mod']
	FP_WinSizes = fire_pixel['FP_WinSize']
	FP_land_water_masks = fire_pixel['FP_land_water_mask']
	
	land_water_suffix = ['_OCEAN', '_LAND']
	FP_data_types = [ land_water_suffix[idx] for idx in FP_land_water_masks.astype(int)]
	temp_holder = {}

	for FP_line, FP_sample, FP_WinSize, data_type in zip(FP_lines, FP_samples, FP_WinSizes, FP_data_types):

		bg_rad_M13  = get_local_rad13(rad_BG, FP_line, FP_sample, FP_WinSize, data_type, namelist, modData)

		for key in bg_rad_M13.keys():
			if key not in temp_holder.keys():
				temp_holder[key] = []
				temp_holder[key].append( bg_rad_M13[key] )
			else:
				temp_holder[key].append( bg_rad_M13[key] )

		
	for key in temp_holder.keys():
		fire_pixel[key] = np.array(temp_holder[key])

	fire_pixel['Solar_Zenith'] = modData['solar_zenith'][FP_lines, FP_samples]

# 	print(' - In FILDA:')
# 	print(fire_pixel.keys())
# 	
# 	idx = np.where( fire_pixel['DNB_observations'] == fire_pixel['DNB_observations'])
# 	
# 	
# 	print(fire_pixel['FP_MDNB_Mean'][idx])
# 	print(fire_pixel['DNB_observations'][idx])
# 	print(fire_pixel['solar_zenith'][idx])
	
	return fire_pixel

#-----------------------------------------------------------------------

def get_fire_paras(fire_pixel, MODAREALUT, namelist):

	'''
	get_fire_paras calculate the fire parameters for each detected fire pixel

	-----------
	Parameters:
	fire_pixel	: dictionary for holding information of the detected fire pixel
	MODAREALUT	: M-band aera look up table
	
	-----------
	Returns:
	fire_pixel	:  dictionary for holding information of the detected fire pixel
	
	'''


	lut 	   = Dataset(MODAREALUT, 'r')	
	area       = lut.variables['MBand_Pixel_Areas'][:]
	lut.close()

	
	# get the unique M band fire pixel
	FP_line_mod   = fire_pixel['FP_line_mod']
	FP_sample_mod = fire_pixel['FP_sample_mod']
	
	num_fire      = []
	for line_mod, sample_mod in zip(FP_line_mod, FP_sample_mod):
	
		add_idx = abs(line_mod - fire_pixel['FP_line_mod']) + abs(sample_mod - fire_pixel['FP_sample_mod'])
		
		idx     = np.where(add_idx == 0)
		
		num_fire.append(np.shape(idx)[1])

	num_fire    = np.array(num_fire)
	FP_Num_Fire = num_fire
	FP_Area_MOD = area[FP_line_mod, FP_sample_mod]
	
	print(f' - FILDA: Calculating fire paramters...' )
	# main part of calculating the fire parameter...
	c     = 2.88 * 10**-9
	sigma = 5.6704*10**-8
	frp   = FP_Area_MOD * sigma * (fire_pixel['FP_M13_Rad'] - fire_pixel['FP_M13_Rad_Mean']) / c * 10**-6

	# sep make a mistake here...
	ve    = ( fire_pixel['DNB_observations'] * 1e-9 * FP_Area_MOD  * np.pi)*(10**4) / 10**6
	
	idx_twilight = np.where( (fire_pixel['Solar_Zenith']>namelist['twilight_ang_min']) & \
						     (fire_pixel['Solar_Zenith']<namelist['twilight_ang_max']) )
	
	ve[idx_twilight] = ( (fire_pixel['DNB_observations'][idx_twilight] -  fire_pixel['FP_DNB_Rad_Mean'][idx_twilight]  )  \
	                      * 1e-9 * FP_Area_MOD[idx_twilight]  * np.pi)*(10**4) / 10**6
	vef   = ve/frp
	
	idx      = np.where(vef!=vef)
	vef[idx] = -999
	
	idx   	 = np.where((frp <0) | (vef<0)  ) #| 
	vef[idx] = 999
	 
	mce   	 = np.log(2*vef)*0.017 + 1
	mce[idx] = np.nan
	vef[idx] = np.nan
	ve[idx]  = np.nan
	frp[idx] = np.nan

	fire_pixel['FP_MCE']		= mce
	fire_pixel['FP_VEF']		= vef
	fire_pixel['FP_Power']		= frp
	fire_pixel['FP_VE']			= ve
	fire_pixel['FP_Area_mod']	= FP_Area_MOD
	fire_pixel['FP_Num_Fire']	= FP_Num_Fire

	
	# remove the pixel if there is no FRP retrieved...
	valid_idx = np.where( (fire_pixel['FP_Power_QA'] != -1) & \
	                      (fire_pixel['FP_Power'] == fire_pixel['FP_Power']))
	
	for key in fire_pixel.keys():
		fire_pixel[key] = fire_pixel[key][valid_idx]
	
	# here we do the SAA check
	
	fire_pixel = check_saa_2(fire_pixel)
	
	return fire_pixel
	
#-----------------------------------------------------------------------
def check_saa_2(fire_pixel):

	delta_M11 = FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M11_Rad'], FILDA_BT.lamda['lamda_' + 'M11']) - \
				FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M11_Rad_Mean'] + 3*fire_pixel['FP_M11_Rad_STD'], FILDA_BT.lamda['lamda_' + 'M11'])

	delta_M12 = FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M12_Rad'], FILDA_BT.lamda['lamda_' + 'M12']) - \
				FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M12_Rad_Mean'] + 3*fire_pixel['FP_M12_Rad_STD'], FILDA_BT.lamda['lamda_' + 'M12'])

	delta_M13 = FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M13_Rad'], FILDA_BT.lamda['lamda_' + 'M13']) - \
				FILDA_BT.cal_brightness_temperature(fire_pixel['FP_M13_Rad_Mean'] + 3*fire_pixel['FP_M13_Rad_MAD'], FILDA_BT.lamda['lamda_' + 'M13'])	

	idx = np.where( ((delta_M11<10) | (delta_M13<2.5)) & (fire_pixel['SAA_flag'] == 1))
	print(f' - Check SAA: Find {np.shape(idx)[1]} SAA affected pixel' )
	SAA_flag = np.zeros_like(fire_pixel['SAA_flag'])
	SAA_flag[idx] = 1
	fire_pixel['FP_SAA_flag'] = SAA_flag
	
	return fire_pixel

#-----------------------------------------------------------------------
def get_surface_type_sinu(fire_img, fire_mod, namelist, time, **kwargs):
	import numpy as np
	import glob
	'''
	get_surface_type prepares the land surface type for each fire pixels

	-----------
	Parameters:

	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	namelist 	: namelist dictionary
	time		: time string to determine the yeas for the land surface type

	-----------
	Returns: 
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution

	'''

	# path of the rad file
	numCeil	  = kwargs.get('numCeil', 2400)
	pixel_size = kwargs.get('pixel_size', 926.62543305)


	SURF_DIR = namelist['SURF_DIR']
	resol = pixel_size/2.0
	year = time.Y
  

	if float(time.Y) >= namelist['Default_year_SURF']:
		year = str(int(namelist['Default_year_SURF']))
		print(f' - Use the default land surface database [MODIS {namelist["Default_year_SURF"]}].')
	else:
		year = time.Y

	if namelist['COPY_STATIC'] == 1:
		SURF_DIR  = SURF_DIR
	else:
		SURF_DIR  = SURF_DIR + year + '/'
	
	# MZ, modify the logic for find the max, min tiles in sinusoidal projection...
	cord = [ fire_img['FP_Latitude'], fire_img['FP_Longitude']]
	tiles, hidMax, hidMin, vidMax, vidMin = FILDA_Time_Cord.get_tile_sinusoidal_3(cord)

	hids = []
	vids = []
	for tile in tiles:
		hids.append( int(float(tile[1:3])) )
		vids.append( int(float(tile[4:])) )
	hids = list(set(hids))
	vids = list(set(vids))    

	#
	GridDim = ((len(vids)) * numCeil, (len(hids)) * numCeil)
	land_type = np.full(GridDim, np.nan)
	mesh_X = np.full(GridDim, np.nan)
	mesh_Y = np.full(GridDim, np.nan)   

	for tile in tiles:
		hid = int(float(tile[1:3]))
		vid = int(float(tile[4:]))
	
		# calculate the coordinates under sinusodal coordinates..
		# they are linear space grids
		x, y = FILDA_Time_Cord.cal_sinu_xy(tile, numCeil)

		hIdx = hid - hidMin 
		vIdx = vid - vidMin

		mesh_X[vIdx * numCeil : (vIdx + 1) * numCeil, \
				hIdx * numCeil : (hIdx + 1) * numCeil, ] = x

		mesh_Y[vIdx * numCeil : (vIdx + 1) * numCeil, \
				hIdx * numCeil : (hIdx + 1) * numCeil, ] = y
	
		# Then try to find the land surface data...
		# Improvement needs to be make here for after MODIS era...
		filename = glob.glob(SURF_DIR + '*' + tile + '*.hdf')    
		if len(filename) == 0:
			continue

		# MZ May-10-2023, define new function to read the land surface type
		# using pyHDF API
		nc_data = FILDA_IO.read_MCD12Q1( filename[0], ['LC_Type1', 'LW'])
	
# 		nc_data = FILDA_IO.read_nc( filename[0], ['LC_Type1', 'LW'])
	
		land_tile = nc_data['LC_Type1']
		water_flag = nc_data['LW']
	
		land_tile[water_flag==1]=17

		land_type[vIdx * numCeil : (vIdx + 1) * numCeil, \
				  hIdx * numCeil : (hIdx + 1) * numCeil, ] = land_tile

	# calculate the boundary of the grid
	x_min = np.nanmin(mesh_X)
	y_max= np.nanmax(mesh_Y)

	# for I band...
	# MZ April. 2023, modify to aviod bugs
	x_s, y_s = FILDA_Time_Cord.geog_to_sinu([fire_img['FP_Latitude'], fire_img['FP_Longitude']])
	x_id = (x_s - x_min + resol/2.)//resol
	y_id = (y_max - y_s + resol/2.)//resol


	x_id = x_id.astype(int)
	y_id = y_id.astype(int)

	fire_img['FP_Land_Type'] = land_type[y_id, x_id]


	# for M band...
	x_s, y_s = FILDA_Time_Cord.geog_to_sinu([fire_mod['FP_Latitude'], fire_mod['FP_Longitude']])
	x_id = (x_s - x_min + resol/2.)//resol + 1
	y_id = (y_max - y_s + resol/2.)//resol + 1

	x_id = x_id.astype(int)
	y_id = y_id.astype(int)

	fire_mod['FP_Land_Type'] = land_type[y_id, x_id]
	
	return fire_img, fire_mod


#-----------------------------------------------------------------------
def get_surface_type(fire_img, fire_mod, namelist, time, **kwargs):
	'''
	get_surface_type prepares the land surface type for each fire pixels
	
	-----------
	Parameters:
	
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	namelist 	: namelist dictionary
	time		: time string to determine the yeas for the land surface type
	
	-----------
	Returns: 
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	
	'''

	# path of the rad file
	numCeil	  = kwargs.get('numCeil', 2400)
	tile_size = kwargs.get('tile_size', 10.)
	SURF_DIR  = namelist['SURF_DIR']
	
	
	year = time.Y
	
	if year not in ['2019', '2018', '2017', '2016', '2015']:
		year = '2018'
		print(' - Use the default land surface database [MODIS 2018].')	
	
	
	SURF_DIR = SURF_DIR + year + '/'
	
	# get the boundary of the geo location...
	north = np.nanmax( (np.nanmax( fire_img['FP_Latitude']  ), np.nanmax( fire_mod['FP_Latitude']  )) )
	south = np.nanmin( (np.nanmin( fire_img['FP_Latitude']  ), np.nanmin( fire_mod['FP_Latitude']  )) )
	west  = np.nanmin( (np.nanmin( fire_img['FP_Longitude'] ), np.nanmin( fire_mod['FP_Longitude'] )) )
	east  = np.nanmax( (np.nanmax( fire_img['FP_Longitude'] ), np.nanmax( fire_mod['FP_Longitude'] )) )
	
	cord  = [north, south, west, east]
	resol = tile_size/numCeil
	# get the tile names for the give boundary
	tiles = FILDA_Time_Cord.get_tiles(cord)
	
	# create the coordinates for the standard grid
	meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = FILDA_Time_Cord.get_cord_PlateCarree(cord, numCeil)	
	
	resol = 10. / numCeil
	
	# create the land type variables...
	land_type = np.full_like(meshLat, np.nan)
	land_type = land_type.astype(int)
	
	for tile in tiles:
		hh = np.int( np.float( tile[1:3] ) )
		vv = np.int( np.float( tile[4:]  ) )

		filename = SURF_DIR +  tile + '.nc'
		if os.path.exists(filename):
			print(' - Reading', filename)
			nc_data    =  FILDA_IO.read_nc( filename, ['land_type'])

			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			for key in nc_data.keys():
				land_type[vIdx * numCeil : (vIdx + 1) * numCeil, \
					 	  hIdx * numCeil : (hIdx + 1) * numCeil, ] = np.flipud(nc_data['land_type'].astype(int))
	
	
	lat_max = np.max(meshLat)
	lat_min = np.min(meshLat)
	
	lon_max = np.max(meshLon)
	lon_min = np.min(meshLon)
	
	# for I band
	latIdx = (lat_max - fire_img['FP_Latitude'])//resol
	latIdx = latIdx.astype(int)
	lonIdx = (fire_img['FP_Longitude'] -lon_min)//resol
	lonIdx = lonIdx.astype(int)

	fire_img['FP_Land_Type'] = land_type[latIdx, lonIdx]

	# for M band
	latIdx = (lat_max - fire_mod['FP_Latitude'])//resol
	latIdx = latIdx.astype(int)
	lonIdx = (fire_mod['FP_Longitude'] -lon_min)//resol
	lonIdx = lonIdx.astype(int)

	fire_mod['FP_Land_Type'] = land_type[latIdx, lonIdx]

	return fire_img, fire_mod

#-----------------------------------------------------------------------
def add_aux_infor(fire_pixel, cdt_fire, namelist, time, modData, imgData, **kwargs):
	
	'''
	add_aux_infor integrates all auxiliary information for the fire pixel
	
	-----------
	Parameters:
	fire_pixel	: dictionary for holding information of the detected fire pixel
	cdt_fire  	: dictionary for holding the fire candidates
	namelist 	: namelist dictionary
	time		: time string to determine the yeas for the land surface type
	modData		: dictionary for holding the M-band data
	imgData		: dictionary for holding the I-band data
	
	-----------
	Returns: 
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	'''
	
	
	flag_land_surface = kwargs.get('flag_land_surface', True)
	flag_gas_flaring  = kwargs.get('flag_gas_flaring', True)
	flag_peatland     = kwargs.get('flag_peatland', True)
	
	# first need to split to I-band and M-BAND, also regularize the name 
	# of the variable...
	fire_img = {}
	fire_mod = {}
	commonVar = ['DNB_observations', 'FP_MCE', 'FP_VEF', 'FP_Status', 'FP_Num_Fire',
				 'FP_I04_Mean', 'FP_I05_Mean', 'FP_BTD_Mean', 'FP_WinSize',
				 'FP_M13_Rad',  'FP_M13_Rad_Mean', 'FP_M13_Rad_MAD', 
				 'FP_M13_Rad_Num', 'FP_M13_WinSize', 'FP_Power_QA',
				 'FP_M07_Rad', 'FP_M07_Rad_Mean', 'FP_M07_Rad_Num',
				 'FP_M08_Rad', 'FP_M08_Rad_Mean', 'FP_M08_Rad_Num',
				 'FP_M10_Rad', 'FP_M10_Rad_Mean', 'FP_M10_Rad_Num',
				 'FP_M11_Rad', 'FP_M11_Rad_Mean', 'FP_M11_Rad_Num',
				 'FP_M12_Rad', 'FP_M12_Rad_Mean', 'FP_M12_Rad_Num',
				 'FP_M14_Rad', 'FP_M14_Rad_Mean', 'FP_M14_Rad_Num',
				 'FP_M15_Rad', 'FP_M15_Rad_Mean', 'FP_M15_Rad_Num',
				 'FP_M16_Rad', 'FP_M16_Rad_Mean', 'FP_M16_Rad_Num',
				 'FP_I04_Rad', 'FP_I04_Rad_Mean', 'FP_I04_Rad_Num', 
				 'FP_I05_Rad', 'FP_I05_Rad_Mean', 'FP_I05_Rad_Num', 
				 'FP_BG_Temp', 'FP_Fire_Temp', 'FP_Fire_Frac', 'FP_Opt_Status']

	fire_img = {}
	for var in commonVar:
		fire_img[var] = fire_pixel[var]
	fire_img['FP_DNB_POS']  = fire_pixel['FP_posDNB']		
	fire_img['FP_Power'] 	= fire_pixel['FP_Power'] / fire_pixel['FP_Num_Fire']
	fire_img['FP_VE']  		= fire_pixel['FP_VE']    / fire_pixel['FP_Num_Fire']
	fire_img['FP_Area']		= fire_pixel['FP_Area_mod']/4.0
	
	fire_img['FP_Line']		= fire_pixel['FP_line_img']
	fire_img['FP_Sample']	= fire_pixel['FP_sample_img']
	fire_img['FP_Latitude']	= fire_pixel['FP_latitude_img']
	fire_img['FP_Longitude']= fire_pixel['FP_longitude_img']
	
	fire_img['FP_IMG_BTD']	= fire_pixel['BTD_IMG']
	fire_img['FP_I04_BT']	= fire_pixel['BTI04']
	fire_img['FP_I05_BT']	= fire_pixel['BTI05']
	fire_img['FP_CM']		= fire_pixel['CM_IMG']
	fire_img['FP_Status']	= fire_pixel['FP_Status']

	fire_img['FP_I04_MAD']	= fire_pixel['FP_I04_MAD']
	fire_img['FP_I05_MAD']	= fire_pixel['FP_I05_MAD']
	fire_img['FP_BTD_MAD']	= fire_pixel['FP_BTD_MAD']
	fire_img['FP_Bowtie']   = 1.0 - get_bowtie(fire_img, 'IMG', namelist)
	fire_img['FP_land_water_mask']	= fire_pixel['FP_land_water_mask']
	fire_img['FP_SAA_flag']			= fire_pixel['FP_SAA_flag']
	
	
	fire_img['Sensor_Zenith']  = imgData['sensor_zenith'][fire_pixel['FP_line_img'], fire_pixel['FP_sample_img']]
	fire_img['Sensor_Azimuth'] = imgData['sensor_azimuth'][fire_pixel['FP_line_img'], fire_pixel['FP_sample_img']]
	fire_img['Fire_mask']      = gen_fire_img_mask(namelist, imgData, fire_img, cdt_fire)

	fire_img['Algorithm_QA']   = gen_algorithm_QA(namelist, imgData, modData, fire_img, cdt_fire)


	fire_img['FP_confidence'] = fire_img['Fire_mask'][fire_pixel['FP_line_img'], fire_pixel['FP_sample_img']]

	fire_img['Solar_Zenith'] = fire_pixel['Solar_Zenith']

	#-------------------------------------------------------------------
	fire_mod = {}
	for var in commonVar:
		fire_mod[var] = fire_pixel[var]
	fire_mod['FP_DNB_POS']  = fire_pixel['FP_posDNB']		
	fire_mod['FP_Power'] 	= fire_pixel['FP_Power']
	fire_mod['FP_VE']  		= fire_pixel['FP_VE']
	fire_mod['FP_Area']		= fire_pixel['FP_Area_mod']
	
	fire_mod['FP_Line']		= fire_pixel['FP_line_mod']
	fire_mod['FP_Sample']	= fire_pixel['FP_sample_mod']
	fire_mod['FP_Latitude']	= fire_pixel['FP_latitude_mod']
	fire_mod['FP_Longitude']= fire_pixel['FP_longitude_mod']
	fire_mod['FP_CM']	    = fire_pixel['CM_MOD']
	fire_mod['FP_Bowtie']   = 1.0 - get_bowtie(fire_mod, 'MOD', namelist)
	
	fire_mod['Solar_Zenith'] = fire_pixel['Solar_Zenith']

	# get the unique fire index...
	fire_pix_mod = np.concatenate( ( np.expand_dims(fire_pixel['FP_line_mod'],   axis = 1), \
									 np.expand_dims(fire_pixel['FP_sample_mod'], axis = 1) ), axis = 1)
	_, unique_idx_mod =  np.unique(fire_pix_mod, return_index=True, axis = 0)
	
	# at last, remove the redundant M band pixel
	for key in fire_mod.keys():
		fire_mod[key] = fire_mod[key][unique_idx_mod]

	# MZ, Mar 26 2024, change M-band fire mask to M-band resolution
	# fire_mod['Fire_mask']   = fire_img['Fire_mask']
	fire_mod['Fire_mask']     = gen_fire_mod_mask(fire_img['Fire_mask'])
	fire_mod['FP_confidence'] = fire_mod['Fire_mask'][fire_mod['FP_Line'], fire_mod['FP_Sample']]	
	fire_mod['Algorithm_QA']  = fire_img['Algorithm_QA']

	print(' - FILDA: Detect', len(fire_img['FP_Line']), 'fires on I band')
	print(' - FILDA: Detect', len(fire_mod['FP_Line']), 'fires on M band')
	
	
	# get the land surface type
	# MZ 10/02/2022, updates to read the MCD12Q1
	# fire_img, fire_mod = get_surface_type(fire_img, fire_mod, namelist, time, **kwargs)
	if flag_land_surface:
		fire_img, fire_mod = get_surface_type_sinu(fire_img, fire_mod, namelist, time)
	
	# get the gas flaring type
	if flag_gas_flaring:
		fire_img, fire_mod = get_gasflaring(fire_img, fire_mod, namelist)
	
	# get the peatland information
	if flag_peatland:
		fire_img, fire_mod = get_peatland(fire_img, fire_mod, namelist)


	fire_img = gen_AdjWater(fire_img, imgData)
	fire_img = gen_AdjCloud(fire_img, imgData)

	fire_mod = gen_AdjWater(fire_mod, modData)
	fire_mod = gen_AdjCloud(fire_mod, modData)

	# then add I band information to M band detection
	mod_line   = fire_mod['FP_Line']
	mod_sample = fire_mod['FP_Sample']


	# MZ 04/23/2024, Add M-band SAA flag
	saa_temp = np.zeros_like(modData['latitude']).astype(int)
	saa_temp[fire_pixel['FP_line_mod'], fire_pixel['FP_sample_mod']] = fire_pixel['FP_SAA_flag']
	fire_mod['FP_SAA_flag'] = saa_temp[mod_line, mod_sample]	

	fire_mod['FP_T04_1'] = imgData['BTI04'][mod_line*2,   mod_sample*2]
	fire_mod['FP_T04_2'] = imgData['BTI04'][mod_line*2+1, mod_sample*2]
	fire_mod['FP_T04_3'] = imgData['BTI04'][mod_line*2,   mod_sample*2+1]
	fire_mod['FP_T04_4'] = imgData['BTI04'][mod_line*2+1, mod_sample*2+1]
	
	fire_mod['FP_T05_1'] = imgData['BTI05'][mod_line*2,   mod_sample*2]
	fire_mod['FP_T05_2'] = imgData['BTI05'][mod_line*2+1, mod_sample*2]
	fire_mod['FP_T05_3'] = imgData['BTI05'][mod_line*2,   mod_sample*2+1]
	fire_mod['FP_T05_4'] = imgData['BTI05'][mod_line*2+1, mod_sample*2+1]
	
	fire_mod['Sensor_Zenith'] = modData['sensor_zenith'][mod_line, mod_sample]
	fire_mod['Sensor_Azimuth'] = modData['sensor_azimuth'][mod_line, mod_sample]
	
	
	return fire_img, fire_mod


#-----------------------------------------------------------------------
def get_gasflaring(fire_img, fire_mod, namelist):
	
	'''
	get_gasflaring prepares the gas flaring flag for the detected fires
	
	-----------
	Parameters:
	
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	namelist	: string, years of the of climatologies 	

	-----------
	Returns: 
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 

	'''

	# path of the rad file
	GAS_DIR  = namelist['GASFLARE_DIR']
	
	# read gas flaring database
# 	gas_file  = 'gas_flaring_2012_2020.nc'
	gas_file  = glob.glob(GAS_DIR + '*.nc')[0]
	print(f' - get_gasflaring: {GAS_DIR}')
	
	gas_ncid = Dataset(gas_file, 'r')
	gas_data = {}
	gas_data['latitude']     = gas_ncid['latitude'][:]
	gas_data['longitude']    = gas_ncid['longitude'][:]
	gas_data['gas_flaring']  = gas_ncid['gas_flaring'][:]
	gas_data['Coordinate']   = gas_ncid.Coordinate
	gas_data['lat_resol']    = float(gas_ncid.Latitude_resolution)
	gas_data['lon_resol']    = float(gas_ncid.Longitude_resolution)
	gas_ncid.close()

	N, S, W, E, _, _, _, _,_ = gas_data['Coordinate'].split(' ')
	S = float(S)
	W = float(W)
	
	gas_data['gas_flaring'][np.where(gas_data['gas_flaring'] != gas_data['gas_flaring'])] = 0

	# find the MOD
	idx_lon = np.round((fire_mod['FP_Longitude'] - W )/gas_data['lon_resol']).astype(int)
	idx_lat = np.round((fire_mod['FP_Latitude'] - S)/gas_data['lat_resol']).astype(int)
	gas_flaring = gas_data['gas_flaring'][idx_lat, idx_lon]
	fire_mod['FP_Gas_Flaring'] = gas_flaring


	# find the IMG
	idx_lon = np.round((fire_img['FP_Longitude'] - W )/gas_data['lon_resol']).astype(int)
	idx_lat = np.round((fire_img['FP_Latitude'] - S)/gas_data['lat_resol']).astype(int)
	gas_flaring = gas_data['gas_flaring'][idx_lat, idx_lon]
	fire_img['FP_Gas_Flaring'] = gas_flaring
	

	return fire_img, fire_mod

#-----------------------------------------------------------------------
def get_peatland(fire_img, fire_mod, namelist, **kwargs):

	'''
	get_peatland prepares the peatland flag for the detected fires
	
	-----------
	Parameters:
	
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 
				  in M-band resolution
	namelist	: string, years of the of climatologies 	

	-----------
	Returns: 
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	fire_mod  	: dictionary for holding information of the detected fire pixel 

	'''
	
	### import numpy as np
	### import os
	
	numCeil	  = kwargs.get('numCeil', 2000)
	tile_size = kwargs.get('tile_size', 10.)


	PEATLAND_DIR  = namelist['PEATLAND_DIR']
	# MZ, the prefix of the peatland now can be specified by the namelist...
	PEATLAND_PREFIX = namelist['PEAT_prefix']
	if PEATLAND_PREFIX == 'N/A':
		print(f' - get_peatland: No prefix for peatland is specified...')
		PEATLAND_PREFIX = ''
	
	# get the boundary of the geo location...
	north = np.nanmax( (np.nanmax( fire_img['FP_Latitude']  ), np.nanmax( fire_mod['FP_Latitude']  )) )
	south = np.nanmin( (np.nanmin( fire_img['FP_Latitude']  ), np.nanmin( fire_mod['FP_Latitude']  )) )
	west  = np.nanmin( (np.nanmin( fire_img['FP_Longitude'] ), np.nanmin( fire_mod['FP_Longitude'] )) )
	east  = np.nanmax( (np.nanmax( fire_img['FP_Longitude'] ), np.nanmax( fire_mod['FP_Longitude'] )) )
	
	cord  = [north, south, west, east]

	resol = tile_size/numCeil
	# get the tile names for the give boundary
	tiles = FILDA_Time_Cord.get_tiles(cord)
	
	# create the coordinates for the standard grid
	meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = FILDA_Time_Cord.get_cord_PlateCarree(cord, numCeil)	
	
	resol = 10. / numCeil

	# create the land type variables...
	peatland  = np.full_like(meshLat, np.nan)
	peatland  = peatland.astype(int)
	peat_fac  = np.full_like(meshLat, np.nan)
	
	for tile in tiles:
		hh = np.int( np.float( tile[1:3] ) )
		vv = np.int( np.float( tile[4:]  ) )

		filename = glob.glob(PEATLAND_DIR + '*' + PEATLAND_PREFIX + '*' + tile + '.nc')
		
		if len(filename)>0:
			print(' - Reading', filename[0])
			nc_data    =  FILDA_IO.read_nc( filename[0], ['area', 'peat_frac'])

			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			for key in nc_data.keys():
				peatland[vIdx * numCeil : (vIdx + 1) * numCeil, \
					 	 hIdx * numCeil : (hIdx + 1) * numCeil, ] = np.flipud(nc_data['area'].astype(int)) /  np.flipud(nc_data['area'].astype(int))

				peat_fac[vIdx * numCeil : (vIdx + 1) * numCeil, \
					 	 hIdx * numCeil : (hIdx + 1) * numCeil, ] = np.flipud(nc_data['peat_frac'])

	peatland[np.where(peatland<0)] = 0
	peat_fac[np.where(peat_fac<0)] = -999

	lat_max = np.max(meshLat)
	lat_min = np.min(meshLat)
	
	lon_max = np.max(meshLon)
	lon_min = np.min(meshLon)

	# for I band
	latIdx = (lat_max - fire_img['FP_Latitude'])//resol
	latIdx = latIdx.astype(int)
	lonIdx = (fire_img['FP_Longitude'] -lon_min)//resol
	lonIdx = lonIdx.astype(int)
	
	# for I band
	latIdx = (lat_max - fire_img['FP_Latitude'])//resol
	latIdx = latIdx.astype(int)
	lonIdx = (fire_img['FP_Longitude'] -lon_min)//resol
	lonIdx = lonIdx.astype(int)

	fire_img['FP_Peatland'] = peatland[latIdx, lonIdx]
	fire_img['FP_Peatfrac'] = peat_fac[latIdx, lonIdx]
	
	# for M band
	latIdx = (lat_max - fire_mod['FP_Latitude'])//resol
	latIdx = latIdx.astype(int)
	lonIdx = (fire_mod['FP_Longitude'] -lon_min)//resol
	lonIdx = lonIdx.astype(int)

	fire_mod['FP_Peatland'] = peatland[latIdx, lonIdx]
	fire_mod['FP_Peatfrac'] = peat_fac[latIdx, lonIdx]

	return fire_img, fire_mod
	

#-----------------------------------------------------------------------
def get_bowtie(fire_pixel, resol, namelist):

	'''
	get_bowtie prepares the bowtie flag for the detected fires
	
	-----------
	Parameters:
	
	fire_pixel  : dictionary for holding information of the detected fire pixel
	resol  		: string to specify the M-band or I-band
	namelist	: string, years of the of climatologies 	

	-----------
	Returns: 
	overlap		: overlap ratio for each fire pixel

	'''

	sat = namelist['platform']
	lut_dir   = namelist['LUT_DIR']
	lut_name  = sat + '_DNB2' + resol + '_Resampling_Lookup_Table.nc'
	
	lut_data  = Dataset(lut_dir + lut_name, 'r')	
	overlap   = lut_data.variables[resol[0] + 'Band_Overlapped'][:]
	overlap   = overlap[fire_pixel['FP_Line'], fire_pixel['FP_Sample']].flatten()
	
	lut_data.close()

	return overlap

#-----------------------------------------------------------------------
def gen_AdjWater(fire_data, lev1bData):

	'''
	gen_AdjWater prepares the water adjacent flag for the detected fires
	
	-----------
	Parameters:
	
	fire_data   : dictionary for holding information of the detected fire pixel
	lev1bData  	: dictionary for holding M/I-band level1b Data
	
	-----------
	Returns: 
	fire_data	: dictionary for holding information of the detected fire pixel

	'''

	land_water_mask = copy.deepcopy(lev1bData['land_water_mask'])
	
	idx = np.where(land_water_mask!=1)  
	land_water_mask[idx] = 0
	land_water_mask = 1 - land_water_mask
	
	# upper left, upper mid, upper right,  mid left, center, mid right, lower left, lower center, lower right
	mask    = land_water_mask[0:-2, 0:-2] + \
			  land_water_mask[0:-2, 1:-1] + \
			  land_water_mask[0:-2, 2:  ] + \
			  land_water_mask[1:-1, 0:-2] + \
			  land_water_mask[1:-1, 2:  ] + \
			  land_water_mask[2: , 0:-2]  + \
			  land_water_mask[2: , 1:-1]  + \
			  land_water_mask[2: , 2:  ]
	
	land_water_mask[1:-1, 1:-1] = land_water_mask[1:-1, 1:-1] + mask
	
	fire_data['FP_AdjWater'] = land_water_mask[fire_data['FP_Line'], fire_data['FP_Sample']]
	return fire_data


#-----------------------------------------------------------------------
def gen_AdjCloud(fire_data, lev1bData):
	
	'''
	gen_AdjCloud prepares the cloud adjacent flag for the detected fires
	
	-----------
	Parameters:
	
	fire_data   : dictionary for holding information of the detected fire pixel
	lev1bData  	: dictionary for holding M/I-band level1b Data
	
	-----------
	Returns: 
	fire_data	: dictionary for holding information of the detected fire pixel

	'''
	CM = copy.deepcopy(lev1bData['CM'])
	
	idx     = np.where(CM!=1)  
	CM[idx] = 0
	CM      = 1 - CM
	# upper left, upper mid, upper right,  mid left, center, mid right, lower left, lower center, lower right
	mask    = CM[0:-2, 0:-2] + \
			  CM[0:-2, 1:-1] + \
			  CM[0:-2, 2:  ] + \
			  CM[1:-1, 0:-2] + \
			  CM[1:-1, 2:  ] + \
			  CM[2:  , 0:-2] + \
			  CM[2:  , 1:-1] + \
			  CM[2:  , 2:  ]
			  
	CM[1:-1, 1:-1] = CM[1:-1, 1:-1] + mask
	

	fire_data['FP_AdjCloud'] = CM[fire_data['FP_Line'], fire_data['FP_Sample']]
	return fire_data

#-----------------------------------------------------------------------
def gen_fire_img_mask(namelist, imgData, fire_img, cdt_fire):

	'''
	gen_fire_img_mask generates the 2D fire mask for the granule data
	
	-----------
	Parameters:
	
	namelist	: namelist dictionary
	imgData		: dictionary for holding the I-band data
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	cdt_fire  	: dictionary for holding the fire candidates
	
	-----------
	Returns: 
	fire_mask	: 2D fire mask

	'''
	
	sat       = namelist['platform']
	lut_dir   = namelist['LUT_DIR']
	lut_name  = sat + '_DNB2IMG_Resampling_Lookup_Table.nc'	
	
	nrow, ncol = np.shape(imgData['latitude'])

	fire_mask = np.full((nrow, ncol), 0)
	
	# get the bowtie information...
	lut_data  = Dataset(lut_dir + lut_name, 'r')
	# MZ, bowtie definition change 
	overlap   = 1.0 - lut_data.variables['IBand_Overlapped'][0:nrow, 0:ncol]      	
	lut_data.close()
	
	# mask water
	idx = np.where(imgData['land_water_mask']!=1)   
	fire_mask[idx] = 3
	
	# mask land
	idx = np.where(imgData['land_water_mask']==1)   
	fire_mask[idx] = 5	
	
	# mask unclassified
	idx = np.where(cdt_fire['QA_FLAG']==-999) 
	idx_x = cdt_fire['FP_line_img'][idx]
	idx_y = cdt_fire['FP_sample_img'][idx]	  
	fire_mask[idx_x, idx_y] = 6	

	# nominal confidence fire pixel 
	idx_x = fire_img['FP_Line']
	idx_y = fire_img['FP_Sample']
	fire_mask[idx_x, idx_y] = 8
	
	# high confidence fire pixel, high confidence fire associated with saturated pixel
	idx   = np.where( fire_img['FP_I04_Mean'] < 0 )
	idx_x = fire_img['FP_Line'][idx]
	idx_y = fire_img['FP_Sample'][idx]
	fire_mask[idx_x, idx_y] = 9

	
	# low confidence fire pixel 
	idx   = np.where((fire_img['FP_land_water_mask'] != 1) & (fire_img['FP_IMG_BTD'] < 15))
	idx_x = fire_img['FP_Line'][idx]
	idx_y = fire_img['FP_Sample'][idx]	
	fire_mask[idx_x, idx_y] = 7
	
	# MZ ARP. 23 2024 follow QA team feedback, restore SAA affect pixel back to 
	# land or ocean
	idx   = np.where((fire_img['FP_SAA_flag'] > 0) & (fire_img['FP_land_water_mask']==1))
	idx_x = fire_img['FP_Line'][idx]
	idx_y = fire_img['FP_Sample'][idx]
	fire_mask[idx_x, idx_y] = 7
	
	
	# MZ follow QA team feedback, restore SAA affect pixel back to 
	# land or ocean
	idx   = np.where((fire_img['FP_SAA_flag'] > 0) & (fire_img['FP_land_water_mask']!=1))
	idx_x = fire_img['FP_Line'][idx]
	idx_y = fire_img['FP_Sample'][idx]
	fire_mask[idx_x, idx_y] = 3

	# mask cloud
	idx = np.where(imgData['CM']<=0)   
	fire_mask[idx] = 4

	# mask bowtie
	idx = np.where(overlap>0.2)
	fire_mask[idx] = 1
	
	# print(1, np.sum(fire_mask>6), fire_mask[441, 799])
	# print('water', np.sum(fire_mask==3), fire_mask[441, 799])
	# print('land', np.sum(fire_mask==5), fire_mask[441, 799])
	
	return fire_mask

#-----------------------------------------------------------------------
def gen_fire_mod_mask(fire_mask_img):

	nrow, ncol = fire_mask_img.shape
	fire_mask_mod = np.full((nrow//2, ncol//2), 0)
	
	# 1 mask land, 5 mask unclassified, 6, mask unclassified
	# 8 nominal confidence fire pixel, 9 high confidence fire pixel
	# 7 low confidence fire pixel, 4 mask cloud
	# 1 mask bowtie
	flags = [3, 5, 6, 8, 9, 7, 4, 1]
	
	for flag in flags:
	# mask water
		idx = np.where(fire_mask_img==flag)
		idx_line, idx_sample = idx[0], idx[1]
		fire_mask_mod[idx_line//2, idx_sample//2] = flag

	return fire_mask_mod
	
#-----------------------------------------------------------------------
def gen_algorithm_QA(namelist, imgData, modData, fire_img, cdt_fire):
	'''
	gen_algorithm_QA generates the 2D algorithm QA
	
	-----------
	Parameters:
	
	namelist	: namelist dictionary
	imgData		: dictionary for holding the I-band data
	modData		: dictionary for holding the M-band data
	fire_img  	: dictionary for holding information of the detected fire pixel 
				  in I-band resolution
	cdt_fire  	: dictionary for holding the fire candidates
	
	-----------
	Returns: 
	QA			: 2D algorithm QA	
	'''

	def reverse_Bits(n, no_of_bits):
		result = 0
		for i in range(no_of_bits):
			result <<= 1
			result |= n & 1
			n >>= 1
		return result

	### from netCDF4 import Dataset
	sat       = namelist['platform']
	lut_dir   = namelist['LUT_DIR']
	lut_name  = sat + '_DNB2IMG_Resampling_Lookup_Table.nc'	
	
	nrow, ncol = np.shape(imgData['latitude'])

	QA   = np.full((nrow, ncol), 0)
	QA   = QA.astype(np.uint32)
	

	# bit 0, Channel I1 quality
	QA = QA 
	
	# bit 1, Channel I2 quality
	QA = QA<<1
	
	# bit 2, Channel I3 quality
	QA = QA<<1

	
	# bit 3, Channel I4 quality	
	validIDX  = np.where( imgData['I04_quality_flags']==0 )
	bit_value = np.full((nrow, ncol), 0)
	bit_value[validIDX] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)
	
	# bit 4, Channel I5 quality	
	validIDX  = np.where( imgData['I05_quality_flags']==0 )
	bit_value = np.full((nrow, ncol), 0)
	bit_value[validIDX] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)
	
	# bit 5, geolocation quality	
	validIDX  = np.where( imgData['quality_flag']==0 )
	bit_value = np.full((nrow, ncol), 0)
	bit_value[validIDX] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 6, Channel M13 quality	
	validIDX  = np.where( modData['M13_quality_flags']==0 )
	FP_lines   = validIDX[0]
	FP_samples = validIDX[1]
	bit_value = np.full((nrow, ncol), 0)
	bit_value[FP_lines, FP_samples] = 1
	bit_value[FP_lines+1, FP_samples] = 1
	bit_value[FP_lines, FP_samples+1] = 1
	bit_value[FP_lines+1, FP_samples+1] = 1
	
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 7, Absolute fire
	bit_value = np.full((nrow, ncol), 0)
	flag      = cdt_fire['FP_abs_img'].astype(int)
	validIDX  = np.where(flag==1)
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1	
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 8, Background high temperature pixel
	bit_value = np.full((nrow, ncol), 0)
	flag      = cdt_fire['FP_bg_img'].astype(int)
	validIDX  = np.where(flag==1)
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1	
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 9, Bright pixel rejection
	QA = QA<<1

	# bit 10, Candidate pixel
	bit_value = np.full((nrow, ncol), 0)
	bit_value[cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 11, Candidate detected through IR
	bit_value = np.full((nrow, ncol), 0)
	flag      = cdt_fire['FP_ir_img'].astype(int)
	validIDX  = np.where(flag==1)
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1	
	QA = QA<<1 | bit_value.astype(np.uint32)


	# bit 12, Candidate detected through visible
	bit_value = np.full((nrow, ncol), 0)
	flag_dnb      = cdt_fire['FP_posDNB']
	flag_abs  = cdt_fire['FP_bg_img'].astype(int)
	flag_bg   = cdt_fire['FP_bg_img'].astype(int)
	flag_ir   = cdt_fire['FP_ir_img']
	validIDX  = np.where((flag_dnb>namelist['thres_DNB']) & (flag_abs != 1) & (flag_bg != 1) & (flag_ir != 1))
	bit_value[cdt_fire['FP_line_img'], cdt_fire['FP_sample_img']] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 13, Scene background
	QA = QA<<1

	# bit 14, Test 1, BTD >= thres_BTD_Dyn
	bit_value = np.full((nrow, ncol), 0)
	flag      = abs(cdt_fire['FP_Status']).astype(int)
	validIDX  = np.where(flag>3)
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 15, Test 2, BTD >= thres_BTD_Abs
	bit_value = np.full((nrow, ncol), 0)
	flag      = abs(cdt_fire['FP_Status']).astype(int)
	validIDX  = np.where( (flag ==7) | (flag ==6)  | (flag ==3) | (flag ==2)  )
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 16, Test 3, I04 >= thres_I04_Dyn
	bit_value = np.full((nrow, ncol), 0)
	flag      = abs(cdt_fire['FP_Status']).astype(int)
	validIDX  = np.where( (flag ==7) | (flag ==6)  | (flag ==3) | (flag ==2)  )
	bit_value[cdt_fire['FP_line_img'][validIDX], cdt_fire['FP_sample_img'][validIDX]] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)
	
	# bit 17, daytime test
	QA = QA<<1

	# bit 18, Saturation...
	validIDX	= np.where( (imgData['BTI05']>325.) | (imgData['BTI04']>355.) | (imgData['BTD_IMG']<0) )
	bit_value = np.full((nrow, ncol), 0)				 
	bit_value[validIDX] = 1								  
	QA = QA<<1 | bit_value.astype(np.uint32)
	
	# bit 19, Glint condition
	QA = QA<<1
	
	# bit 20, South Atlantic...
	QA = QA<<1
	
	# bit 21, fire pixel over water
	bit_value = np.full((nrow, ncol), 0)
	fire_water = imgData['land_water_mask'][fire_img['FP_Line'], fire_img['FP_Sample']]
	validIDX = np.where(fire_water != 1)
	bit_value[fire_img['FP_Line'][validIDX],  fire_img['FP_Sample'][validIDX]] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)

	# bit 22, twilight pixel
	bit_value = np.full((nrow, ncol), 0)
	validIDX = np.where( (imgData['solar_zenith']>namelist['twilight_ang_min']) & \
					     (imgData['solar_zenith']<namelist['twilight_ang_max']) )
	bit_value[validIDX] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)
		
	# bit 23, persistence test temperature
	QA = QA<<1

	# bit 24, persistence Test number 
	QA = QA<<1

	# bit 25, bowtie	
	sat       = namelist['platform']
	lut_dir   = namelist['LUT_DIR']
	lut_name  = sat + '_DNB2IMG_Resampling_Lookup_Table.nc'		
	lut_data  = Dataset(lut_dir + lut_name, 'r')
	# MZ, bowtie definition change 
	overlap   = 1.-lut_data.variables['IBand_Overlapped'][0:nrow, 0:ncol]      	
	idx = np.where(overlap>0.2)
	bit_value = np.full((nrow, ncol), 0)
	bit_value[idx] = 1
	QA = QA<<1 | bit_value.astype(np.uint32)
	
	# bit 26-31
	QA = QA<<4

	QA = reverse_Bits(QA, 32)

	return QA
	
