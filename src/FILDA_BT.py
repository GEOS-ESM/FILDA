'''
Comment in progress
'''

import numpy as np

lamda = dict(lamda_M07 = 0.865,
		 	lamda_M08 = 1.24,
		 	lamda_M10 = 1.61,
		 	lamda_M11 = 2.25,
		 	lamda_M12 = 3.70,
		 	lamda_M13 = 4.05,
		 	lamda_M14 = 8.55,
			lamda_M15 = 10.76,
			lamda_M16 = 12.0125,
	 	 	lamda_I04 = 3.757557,	#3.757494
		 	lamda_I05 = 11.57841)	#11.587094

transmittance =dict(M07 = 0.9851676,  # 0.865,
					M08 = 0.9834966,  # 1.24,
					M10 = 0.9480004,  # 1.61,
					M11 = 0.95913553, # 2.25,
					M12 = 0.89846635, # 3.70,
					M13 = 0.7823044,  # 4.05,
					M14 = 0.68727845, # 8.55,
					M15 = 0.8056259,  # 10.76
					M16 = 0.7417919 )  # 12.0125)


#-----------------------------------------------------------------------
def get_bt(modData, imgData, mod_band = ['M12', 'M13','M15', 'M16'], remove_rad = True):
	
	#import numpy as np
	
	print(' - FILDA: Calculate brightness temperature...')
	for band in mod_band:
		lut_value = np.ma.filled(modData[band+'_LUT'], fill_value = np.nan)
		if np.nansum( lut_value) <=0:
			print(' - FILDA: Calculate brightness temperature using redefined wavelength, band', band)
			modData['BT' + band] = cal_brightness_temperature(modData[band+'_rad'], lamda['lamda_' + band])
	
	# calculate the brightness temperature on I band
	lut_value = np.ma.filled(imgData['I04_LUT'], fill_value = np.nan)
	if np.nansum( lut_value)<=0:
		print(' - FILDA: Calculate brightness temperature using redefined wavelength, band I04')
		imgData['BTI04'] = cal_brightness_temperature(imgData['I04_rad'], lamda['lamda_I04'])	# 03.74 um
		
	lut_value = np.ma.filled(imgData['I05_LUT'], fill_value = np.nan)	
	if np.nansum( lut_value) <=0:
		print(' - FILDA: Calculate brightness temperature using redefined wavelength, band I05')
		imgData['BTI05'] = cal_brightness_temperature(imgData['I05_rad'], lamda['lamda_I05'])	# 11.45 um	

	# MZ, Mar. 20, 2024, add output for GMAO
	temp_array = np.array([imgData['I04_rad'][::2,::2], imgData['I04_rad'][1::2,::2], imgData['I04_rad'][::2,1::2], imgData['I04_rad'][1::2,1::2]])
	modData['I04_rad'] = np.nanmean(temp_array, axis = 0)

	temp_array = np.array([imgData['I05_rad'][::2,::2], imgData['I05_rad'][1::2,::2], imgData['I05_rad'][::2,1::2], imgData['I05_rad'][1::2,1::2]])
	modData['I05_rad'] = np.nanmean(temp_array, axis = 0)

	if remove_rad: 
		print(' - FILDA: Remove radiance data from memory...')
		
		modkeys = list(modData.keys())
		# MZ, Based on GMAO needs, we no longer remove M12 radiance
# 		rad_key = ['M12_rad']
		rad_key = []
		for key in modkeys:
			if key in rad_key:
				modData[key] = None
				modData.pop(key)
		
		modkeys = list(modData.keys())	
		for key in modkeys:
			if ("_refl" in key):
				modData[key] = None
				del modData[key]
				
		imgkeys = list(imgData.keys())
		for key in imgkeys:
			if ("_rad" in key) :
				imgData[key] = None
				del imgData[key]
	
				
	modData['BTD_MOD']   = modData['BTM13'] - modData['BTM15']
	imgData['BTD_IMG']   = imgData['BTI04'] - imgData['BTI05']

	return modData, imgData

#-----------------------------------------------------------------------
def cal_brightness_temperature(rad, lamda):
	'''
	rad: 	unit W/m2/sr/um
	lamda:	um
		
	Reference: https://ncc.nesdis.noaa.gov/data/planck.html
	'''
	
	#import numpy as np
	
	h = 6.62607015e-34		# Plank's const in Js
	k = 1.38064852e-23		# Boltzmanns const in J/K. 
	c = 299792458			# speed of light in m/s
	
	lamda = lamda * 1e-6	# Convert um to m
	rad   = rad * 1e6		# Convert W/sr/m2/um to W/sr/m2/m

	

	temperature = (h * c) / (k * lamda) * \
			      (1 / np.log( 2 * h * c**2 * lamda**(-5) / rad + 1))	
		
	return temperature


#-----------------------------------------------------------------------
def planck_func(lamda,T):
	'''
	rad: 	unit W/m2/sr/um
	lamda:	um
		
	Reference: https://ncc.nesdis.noaa.gov/data/planck.html
	
	
	Watts/meter^2/steradian/micrometer
	
	'''



	#import numpy as np
	
	c1  = 1.191042e8
	c2  = 1.4387752e4

	rad = c1/lamda**5/( np.exp(c2/lamda/T) -1)

	
# 	h = 6.626*10**-34
# 	kb = 1.381*10**-23
# 	c = 2.998*10**8
# 	rad= 2*h*(c**2)/((lamda**5)*(np.exp(h*c/(kb*lamda*T))-1))/1000/1000    
	return rad

	
	
	
	
	
	