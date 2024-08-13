'''
Comment in progress
'''

import numpy as np

def cloud_test(modData, imgData, tempField, namelist):
	
	### import numpy as np
	print(' - FILDA: Do cloud screening...')
	cm_mod 			= cloud_test_mod(modData, tempField, namelist)
	
	modData['CM_2']	= np.full_like(modData['latitude'], 0) 
	
	modData['CM_2'][np.where((cm_mod['cloud_gross']   >=1)   & \
							(cm_mod['cloud_infrared']>=1)   & \
							(cm_mod['cloud_water']   >=1))] = 1


	cm_img 			= cloud_test_img(imgData, namelist)
	imgData['CM']	= cm_img['cm']
	
	# modify the cloud mask, use the I band cloud mask...
	modData['CM'] =  (cm_img['cm'][1::2, 1::2] + cm_img['cm'][0:-1:2, 1::2] + \
					  cm_img['cm'][0:-1:2, 0:-1:2] + cm_img['cm'][1::2, 0:-1:2])/4.
					 
	modData['CM'][np.where(modData['CM']>0)] = 1
	
	return modData, imgData
	
	
#-----------------------------------------------------------------------	
def cloud_test_img(imgData, namelist):
	### import numpy as np
	thres_i4			= namelist['thres_cloud_I04']
	thres_i5			= namelist['thres_cloud_I05']	
	
	idx = np.where(imgData['BTI04'] != imgData['BTI04'] )

	cloudmask_1 = np.ones_like(imgData['BTI04'])
	cloudmask_1[np.where(imgData['BTI04']<thres_i4)] = 0
	cloudmask_1[idx] = 0

	cloudmask_2 = np.ones_like(imgData['BTI05'])
	cloudmask_2[np.where(imgData['BTI05']<thres_i5)] = 0
	cloudmask_2[idx] = 0	
	
	cloudmask_final = np.ones_like(imgData['BTI05'])
	
	cloudmask_final[np.where((imgData['BTI04']<thres_i4)&(imgData['BTI05']<thres_i5))]  = 0
	cloudmask_final[idx]  = 0

	cloudmask = {}
	cloudmask['I04_msk']	= cloudmask_1	
	cloudmask['I05_msk']	= cloudmask_2
	cloudmask['cm']		    = cloudmask_final
	
	return cloudmask
	
#-----------------------------------------------------------------------
def cloud_test_mod(modData, tempField, namelist):
	'''
	VIIRS IR cloud test
	
	parameter:
		mbandData: mband radiance data
		tempField: geos-fp temperature field
		infraredLUT: look up table for infrared data
		threshold_cloud_gross: threshold for gross cloud test
		threshold_cloud_high : threshold for high cloud test
		threshold_cloud_water: threshold for water cloud test
	return:
		cloudmask: dictionary like, contain all the results of cloud test
	
	'''
	### import numpy as np
	#--------------------------
	# some constant
	#--------------------------

	infraredLUT     		= namelist['LUT_DIR'] + 'Infrared.csv'
	threshold_cloud_gross	= namelist['thres_cloud_gross']
	threshold_cloud_high	= namelist['thres_cloud_high']
	threshold_cloud_water	= namelist['thres_cloud_water']
	resol_lat			    = namelist['resol_lat']
	resol_lon			    = namelist['resol_lon']		


	infraredLUT = np.genfromtxt(infraredLUT)

	idx = np.where(modData['BTM12'] != modData['BTM12'] )

	cloudmask_1 = cloud_gross(modData['BTM15'], modData['latitude'], modData['longitude'], \
							  resol_lat, resol_lon, tempField, threshold_cloud_gross)
	cloudmask_1[idx] = 0
									
	cloudmask_2 = cloud_infrared(modData['BTM15'], modData['BTM16'], modData['sensor_zenith'], infraredLUT)
	cloudmask_2[idx] = 0
	
	cloudmask_3 = cloud_high(modData['BTM12'], modData['BTM16'], threshold_cloud_high)
	cloudmask_3[idx] = 0
	
	cloudmask_4 = cloud_water(modData['BTM15'], modData['BTM12'], threshold_cloud_water)
	cloudmask_4[idx] = 0
	
	cloudmask_5 = np.ones_like(modData['latitude'])
	cloudmask_5[np.where(modData['BTM12']<280)] = 0
	cloudmask_5[idx] = 0

	cloudmask_6 = np.ones_like(modData['latitude'])
	cloudmask_6[np.where(modData['BTM16']<260)] = 0
	cloudmask_6[idx] = 0
	
	cloudmask_final = cloudmask_1 * cloudmask_2 * cloudmask_3 * cloudmask_4 * cloudmask_5 * cloudmask_6
	
	cloudmask = {}
	
	cloudmask['cloud_gross']	= cloudmask_1
	cloudmask['cloud_infrared'] = cloudmask_2
	cloudmask['cloud_high']		= cloudmask_3
	cloudmask['cloud_water']	= cloudmask_4
	cloudmask['cloud_bt12']	    = cloudmask_5
	cloudmask['cloud_bt16']	    = cloudmask_6
	cloudmask['cloudmask']		= cloudmask_final
	
	return cloudmask

#----------------------------------------------------------------------------
def cloud_gross(bt, lat, lon, resol_lat, resol_lon, tempField, threshold):
	'''
	cloud test 1
	
	BTM15 test
	
	'''
	### import numpy as np
	# create the surface temperature for the M band observation
	latIdx = (lat + 90)//resol_lat
	lonIdx = (lon + 180)//resol_lon

	invalid_idx = np.where( (lat != lat) | (lon != lon) )
	latIdx[invalid_idx] = 0
	lonIdx[invalid_idx] = 0	
	
	latIdx = latIdx.astype(int)
	lonIdx = lonIdx.astype(int)

	surf_temp = tempField[latIdx,lonIdx]
	surf_temp[latIdx[invalid_idx],lonIdx[invalid_idx]] = np.nan
	
	
	cloudmask = np.full_like(bt, 0)
	temp_dif = surf_temp - bt
	cloudmask[np.where(temp_dif < threshold)] = 1	
	cloudmask[np.where(bt != bt)] = np.nan
	return cloudmask

#----------------------------------------------------------------------------
def cloud_infrared(bt1, bt2, vza, lut):
	
	### import numpy as np
	btd = bt1 - bt2

	# create the threshold based on look up table
	# index of the secant of the view zenith angle
	secant_vza = 1.0 / np.cos(vza/180 * np.pi)
	secant_vza[np.where(secant_vza > 2 )] = 2

	seccantIdx = (secant_vza - 1.0)//0.25
	seccantIdx = seccantIdx.astype(int)
	invalidAng = np.where(seccantIdx<0)
	seccantIdx[invalidAng] = 0

	# index of the secant of the m15 brightness temperature
	bt1[np.where(bt1<190)] = 190
	btIdx = abs(bt1 - 310) // 10
	btIdx = btIdx.astype(int)

	invalidBt = np.where(btIdx<0)
	btIdx[invalidBt] = 0

	threshold =  lut[btIdx, seccantIdx]

	cloudmask 								= np.full_like(btd, 0)
	cloudmask[np.where(btd < threshold)]	= 1
	cloudmask[invalidAng]					= np.nan
	cloudmask[invalidBt]					= np.nan

	return cloudmask

#----------------------------------------------------------------------------
def cloud_high(bt1, bt2, threshold):
	
	### import numpy as np
	
	btd = bt1 - bt2
	
	cloudmask = np.full_like(btd, 0)
	cloudmask[np.where(btd < threshold)] = 1
	cloudmask[np.where(bt1 != bt1)] = np.nan
	
	return cloudmask
	
#----------------------------------------------------------------------------	
def cloud_water(bt1, bt2, threshold):
	
	### import numpy as np
	
	btd = bt1 - bt2
	cloudmask = np.full_like(btd, 0)
	cloudmask[np.where(btd < threshold)] = 1
	cloudmask[np.where(bt1 != bt1)] = np.nan
	return cloudmask
