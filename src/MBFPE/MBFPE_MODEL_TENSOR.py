'''
'''

try:
    # Try importing Theano and PyMC3
    import theano
    import theano.tensor as tensor
    import pymc3 as pm
except ImportError:
    # Fall back to pytensor and pymc if Theano/PyMC3 are not available
    try:
        import pytensor as theano
        import pytensor.tensor as tensor
        import pymc as pm
    except ImportError:
        raise ImportError("Required libraries are not installed.")


import numpy as np
import time
# import aesara
# import aesara.tensor as tensor
# import theano
# import theano.tensor as tensor
# import pymc as pm
# import pytensor.tensor as tensor
# import pymc3 as pm



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def trapezoidal_rule_tensor(y, x):
	'''
	# Basic trapezoidal rule implementation for Aesara tensors
	'''
	dx = x[:,1:, :] - x[:,:-1, :]
	avg_height = (y[:,1:, :] + y[:,:-1, :]) / 2.0
	area = tensor.sum(dx * avg_height, axis = 1)
	return area
	

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def planck_func_3D_tensor(lamda, T):
	"""
	Calculate the spectral radiance using Aesara tensors.
	"""
	c1 = 1.191042e8
	c2 = 1.4387752e4
	
	
	# Reshape lamda to (n, 1, 1) and T to (1, m, p) for broadcasting
	lamda = np.reshape(lamda, lamda.shape + (1,))
	# Convert inputs to Aesara tensors
	T_ae = tensor.reshape(T, newshape=(1, 1, T.shape[0]))
	
	# Calculate spectral radiance
	rad = c1 / lamda**5 / (tensor.exp(c2 / lamda / T_ae) - 1)
	return rad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
def get_band_radiance_BG_tensor(T_b, C, bg_tensor, sensor_config):
	"""
	Function to compute the band radiance of a cool background using 
	tensor operations.
	
	Parameters:
	T_b: Background temperature, in K
	vza: Satellite view zenith angle, in degree
	C: Water vapor mass scaling factor, unitless
	lambdas: wavelength of interest
	tau_wvp: water vapor optical depth, unitless
	tau_other_gas: atmosphere trace optical depth (without WVP), unitless
	emissitivity: Bandwise surface emissitivity, unitless
	rsr: Relative sensor response function, unitless
	
	Return:
	
	band_radiance: band radiance for the given band, temperature, sensor
	characteristics, surface and atmosphere conditions
	
	"""
	
	lambdas = sensor_config.lambdas
	tau_wvp = sensor_config.tau_wvp
	tau_other_gas = sensor_config.tau_other_gas
	rsr = sensor_config.rsr
	
	vza = bg_tensor.vza
	emissitivity = bg_tensor.emit
	
	# Convert inputs totensors
	T = tensor.stack([T_b]).T
	# Calculate the radiance based on a given temperature and wavelength
	rad = planck_func_3D_tensor(lambdas, T)
	
	# some necessary to match the dimension of the array for array 
	# broadcast 
	tau_wvp = np.reshape(tau_wvp, tau_wvp.shape + (1,))
	tau_other_gas = np.reshape(tau_other_gas, tau_other_gas.shape + (1,))
	
	rsr = np.reshape(rsr, rsr.shape + (1,))
	wavelengh = np.reshape(lambdas, lambdas.shape + (1,))

	# Atmoshperic correction
	# 2.96 is the simulated standard mass of the column water vapor content
	# make those variables as tensor
	airmass = tensor.cos( tensor.deg2rad(vza) )
	tau_wvp_tensor = tensor.as_tensor_variable(tau_wvp)/2.96
	tau_tau_other_gas = tensor.as_tensor_variable(tau_other_gas)
	
	# Compute product for integration
	# simple Beer-Lambert Law correction
	product = tensor.exp( - ( C * tau_wvp_tensor + tau_tau_other_gas) / \
			  airmass )  * rsr * rad
	
	# Perform integration using the trapezoidal rule
	band_radiance = trapezoidal_rule_tensor(product, wavelengh) / \
					np.trapz(rsr, wavelengh, axis = 1)
	
	# emissitivity is the broadband and applied to the convoluted band radiance
	emissitivity = tensor.reshape(emissitivity, newshape=(emissitivity.shape[0], 1))

	band_radiance = band_radiance*emissitivity
	# Flatten the band_radiance from (N, 1) to (N,)
	band_radiance = tensor.flatten(band_radiance, 1)
	
	return band_radiance


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
def get_band_biphasic_radiance_tensor(T_s, T_f, F_s, F_f, fire_tensor, sensor_config):
	"""
	Function to compute the band radiance of a biphasic fire model using 
	tensor operations.
	
	Parameters:
	T_s: Smoldering temperature 
	T_f: Flaming temperature 
	F_s: Smoldering fraction 
	F_f: Flaming fraction
	vza: Satellite view zenith angle, in degree
	C: Water vapor mass scaling factor, unitless
	lambdas: wavelength of interest
	tau_wvp: water vapor optical depth, unitless
	tau_other_gas: atmosphere trace optical depth (without WVP), unitless
	emissitivity: Bandwise surface emissitivity, unitless
	rsr: Relative sensor response function, unitless
	
	Return:
	band_radiance: band radiance for the given band, temperature, sensor
	characteristics, surface and atmosphere conditions
	
	"""
	
	lambdas = sensor_config.lambdas
	tau_wvp = sensor_config.tau_wvp
	tau_other_gas = sensor_config.tau_other_gas
	rsr = sensor_config.rsr
	
	vza = fire_tensor.vza
	C = fire_tensor.C
	
	# Convert inputs to Aesara tensors
	T = tensor.stack([T_s, T_f]).T
	
	# Calculate the radiance based on a given temperature and wavelength
	rad = planck_func_3D_tensor(lambdas, T)
	
	tau_wvp = np.reshape(tau_wvp, tau_wvp.shape + (1,))
	tau_other_gas = np.reshape(tau_other_gas, tau_other_gas.shape + (1,))
	
	rsr = np.reshape(rsr, rsr.shape + (1,))
	wavelengh = np.reshape(lambdas, lambdas.shape + (1,))
	
	# Compute product for integration
	airmass = tensor.cos( tensor.deg2rad(vza) )
	tau_wvp_tensor = tensor.as_tensor_variable(tau_wvp)/2.96
	tau_tau_other_gas = tensor.as_tensor_variable(tau_other_gas)
	
	product = tensor.exp( - ( C * tau_wvp_tensor + tau_tau_other_gas) / \
			  airmass )  * rsr * rad
	
	# Perform integration using the trapezoidal rule
	band_radiance = trapezoidal_rule_tensor(product, wavelengh) / np.trapz(rsr, wavelengh, axis = 1)
	
	# Weight by the fire fraction.
	frac = tensor.stack([F_s, F_f])
	
	fire_sig = tensor.sum(band_radiance * frac, axis=1)
	
	return fire_sig

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_band_uniphasic_radiance_tensor(T_mean, F_mean, fire_tensor, sensor_config):
	"""
	Function to compute the band radiance of a uniphasic fire model using 
	tensor operations.
	
	Parameters:
	T_mean
	F_mean
	fire_tensor
	sensor_config
	
	Return:
	band_radiance: band radiance for the given band, temperature, sensor
	characteristics, surface and atmosphere conditions
	
	"""

	lambdas = sensor_config.lambdas
	tau_wvp = sensor_config.tau_wvp
	tau_other_gas = sensor_config.tau_other_gas
	rsr = sensor_config.rsr
	
	vza = fire_tensor.vza
	C = fire_tensor.C
	

	# Convert inputs to Aesara tensors
	T = tensor.stack([T_mean]).T

	# Calculate the radiance based on a given temperature and wavelength
	rad = planck_func_3D_tensor(lambdas, T)
	
	tau_wvp = np.reshape(tau_wvp, tau_wvp.shape + (1,))
	tau_other_gas = np.reshape(tau_other_gas, tau_other_gas.shape + (1,))
	
	rsr = np.reshape(rsr, rsr.shape + (1,))
	wavelengh = np.reshape(lambdas, lambdas.shape + (1,))
	
	# Compute product for integration
	airmass = tensor.cos( tensor.deg2rad(vza) )
	tau_wvp_tensor = tensor.as_tensor_variable(tau_wvp)/2.96
	tau_tau_other_gas = tensor.as_tensor_variable(tau_other_gas)
	
	product = tensor.exp( - ( C * tau_wvp_tensor + tau_tau_other_gas) / \
			  airmass )  * rsr * rad
	
	# Perform integration using the trapezoidal rule
	band_radiance = trapezoidal_rule_tensor(product, wavelengh) / np.trapz(rsr, wavelengh, axis = 1)

	frac = tensor.stack([F_mean])
	
	fire_sig = tensor.sum(band_radiance * frac, axis=1)
	
	return fire_sig

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
def get_bi_frp_tensor(T_s, T_f, F_s, F_f, fire_tensor):
	area = fire_tensor.area
	frp =  area * 5.6704*10**-8 * (T_s**4 * F_s +  T_f**4 * F_f) * 1e-6
	return frp
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
def get_uni_frp_tensor(T_mean, F_mean, fire_tensor):
	area = fire_tensor.area
	frp =  area * 5.6704*10**-8 * T_mean**4 * F_mean * 1e-6
	return frp
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
def get_fire_fraction_tensor(F_s, F_f):
	total_frac =  F_s + F_f
	return total_frac
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
if __name__ == '__main__':
	print('')




