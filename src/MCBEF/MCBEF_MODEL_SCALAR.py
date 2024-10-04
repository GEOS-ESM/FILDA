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
	
	return rad
	
#-----------------------------------------------------------------------
def inv_planck_func(rad, lamda):
	'''
	rad: 	unit W/m2/sr/um
	lamda:	um
		
	Reference: https://ncc.nesdis.noaa.gov/data/planck.html
	'''
	
	import numpy as np
	
	h = 6.62607015e-34		# Plank's const in Js
	k = 1.38064852e-23		# Boltzmanns const in J/K. 
	c = 299792458			# speed of light in m/s
	
	lamda = lamda * 1e-6	# Convert um to m
	rad   = rad * 1e6		# Convert W/sr/m2/um to W/sr/m2/m
	
	
	
	temperature = (h * c) / (k * lamda) * \
				  (1 / np.log( 2 * h * c**2 * lamda**(-5) / rad + 1))	
		
	return temperature

#-----------------------------------------------------------------------
def planck_func_2D(lamda, T):
	"""
	Calculate the spectral radiance.
	
	Parameters:
	lamda (1D array): Wavelengths in micrometers.
	T (2D array): Temperature in Kelvin.
	
	Returns:
	2D array: Spectral radiance in units of W/m^2/sr/um, in order of 
	[lamda, temperature]
	
	"""
	c1 = 1.191042e8
	c2 = 1.4387752e4
	
	# Reshape lamda to (n, 1, 1) and T to (1, m, p) for broadcasting
	lamda = np.reshape(lamda, (-1, 1))
	T = np.reshape(T, (1,) + T.shape)
	
	# Calculate spectral radiance using broadcasting
	rad = c1 / lamda**5 / (np.exp(c2 / lamda / T) - 1)
	
	return rad
	
	
	
	
	
def get_band_radiance(T_s, T_f, T_b, F_s, F_f, vza, C, lambdas, transmittance, rsr):

    # map the variable into array...
    T = np.array([T_s, T_f, T_b])
#     frac = np.array([F_s, F_f, 1-(F_s+F_f)])
    frac = np.array([F_s, F_f, 1])
    rad = planck_func_3D(lambdas, T)

    # expand dimension
    transmittance = np.reshape(transmittance, transmittance.shape + (1,))
    rsr = np.reshape(rsr, rsr.shape + (1,))
    wavelengh = np.reshape(lambdas, lambdas.shape + (1,))
    # get the term in the integration braket
    airmass = np.cos( np.deg2rad(vza) )
    
    C = np.reshape(C, C.shape + (1, 1))
    print(np.shape(transmittance))
    product = np.exp(( np.log(transmittance) * C/ airmass) )  * rsr * rad	
    
    
    # do the integration...
    band_radiance = np.trapz(product, wavelengh, axis = 1)/np.trapz(rsr, wavelengh, axis = 1)

    
    
    # get the signal
    fire_sig = np.sum(band_radiance * frac, axis = 1)

    return fire_sig

	
	
	
	
	
