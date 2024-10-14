'''
This is the main library of MonteCarlo Biphasic Estimation of Fire Properties
(MCBEF)

 * Main developer : Meng Zhou (MZ)
 * Supervisor     : Arlindo da Silva (ADS)

Institution : 
 * Global Modeling and Assimilation Office (GMAO), 
   the National Aeronautics and Space Administration (NASA)
 * Goddard Earth Sciences Technology and Research (GESTAR) II, 
   University of Maryland, Baltimore County


Dependence:

MCBEF dependents on the following python libraries:

#	libraries		version
1	os 				Default  Miscellaneous operating system interfaces
2	sys				Default  System-specific parameters and functions
3	copy			Default  Shallow and deep copy operations
4	datetime		Default  Basic date and time types
9	multiprocessing	Default  Process-based parallelism
5	numpy			1.19.1   https://numpy.org/
6	scipy			1.4.1    https://scipy.org/
7	pandas			1.1.3    https://pandas.pydata.org/		
8	netCDF4			1.5.1.2  https://unidata.github.io/netcdf4-python/

'''

#-----------------------------------------------------------------------
# inport the build libraries used in MCBEF

# import theano as theano
# import theano.tensor as tensor
# import pymc3 as pm

try:
    # Try importing Theano and PyMC3
    import theano
    import theano.tensor as tensor
    import pymc3 as pm
    print(" - MCBEF: Imported theano with pymc3")
    
except ImportError:
    # Fall back to pytensor and pymc if Theano/PyMC3 are not available
    try:
        import pytensor as theano
        import pytensor.tensor as tensor
        import pymc as pm
        print(" - MCBEF: Imported pytensor with pymc")
    except ImportError:
        raise ImportError("Required libraries are not installed.")

pymc_version = pm.__version__
# Split the version string and convert to a tuple of integers
version_info = tuple(map(int, pymc_version.split('.')))

# Set a flag for version
is_pymc3 = version_info[0] == 3



import numpy as np
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Process
from functools import partial
import time
import sys
from scipy import stats
from scipy import interpolate
import copy

from . import MCBEF_IO as IO # module for reading and writing files
from . import MCBEF_SPACE_TIME as ST # module 
from . import MCBEF_SENSOR as SR # module for 
from . import MCBEF_NAMELIST as NL # module for 
from . import MCBEF_MODEL_TENSOR as MT # module for
from . import MCBEF_EVAL as EV # module for 
from .MCBEF_utils import printf

import logging
import warnings

# multiprocessing.set_start_method('fork')

# Suppress informational messages and below (e.g., INFO, DEBUG, WARNING)
# from PyMC3
logging.getLogger('pymc3').setLevel(logging.ERROR)
logging.getLogger('arviz').setLevel(logging.ERROR)

# Suppress runtime warnings (e.g., overflow encountered in exp)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

FLAG_UNIPHASIC    = 1
FLAG_BIPHASIC     = 2
FLAG_DEDGRADATE   = 3
FLAG_BACKGROUND   = 10

FLAG_BOWTIE       = 100
FLAG_MISSING_BG   = 101
FLAG_MISSING_FIRE = 102
FLAG_FAIL_UNIPHASIC = 103
FLAG_FAIL_BIPHASIC  = 104
FILLVALUE         = 254
FILLVALUE_FRP     = np.nan



class data:
    def __init__(self, signal, sigma, vza, area):      
        self.signal = signal
        self.sigma = sigma
        self.vza = vza
        self.area = area
           
    def __len__(self):
        return len(self.signal)
    
    def __getitem__(self, idx):
        return self.signal[idx], self.sigma[idx], self.vza[idx], self.area[idx]

class BackgroundData(data):
	def __init__(self, signal, sigma, vza, area, emit, lst_nig, lst_nig_std, lst_day, lst_day_std):
		# Initialize the parent class with the existing attributes
		super().__init__(signal, sigma, vza, area)
		# Initialize the new attribute specific to this subclass
		self.emit = emit
		self.lst_nig = lst_nig
		self.lst_nig_std = lst_nig_std
		self.lst_day = lst_day
		self.lst_day_std = lst_day_std
	
	def __getitem__(self, idx):
		# Get the basic items from the parent class
		items = super().__getitem__(idx)
		# Add the new 'emit' item and return
		return items + (self.emit[idx], self.lst_nig[idx], 
		                self.lst_nig_std[idx], self.lst_day[idx], 
		                self.lst_day_std[idx])

class FireData(data):
	def __init__(self, signal, sigma, vza, area, FRP, 
					   gasflaring, static_flag, bowtie, raw):
		# Initialize the parent class with the existing attributes
		super().__init__(signal, sigma, vza, area)
		# Initialize the new attribute specific to this subclass
		self.raw = raw
		self.FRP = FRP
		self.gasflaring  = gasflaring
		self.static_flag = static_flag
		self.bowtie = bowtie
	
	def __getitem__(self, idx):
		# Get the basic items from the parent class
		items = super().__getitem__(idx)
		# Add the new 'emit' item and return
		return items + (self.FRP[idx], self.gasflaring[idx], 
		                self.static_flag[idx], self.bowtie[idx])
		

def get_sample_set(filda_dict, nl):
	# process the background signal
	snr_bg = np.array([nl.SNR_bg[band] for band in nl.sel_bg_bands])
	band_lst_bg = ['FP_' + band + '_Rad_Mean' for band in nl.sel_bg_bands]
	sig_bg = filda_dict[band_lst_bg].values
	sigma_bg = sig_bg / np.reshape(snr_bg, (1,) + snr_bg.shape)
	
	band_lst_emit_bg = ['Emis_' + band for band in nl.sel_bg_bands]
	emit_bg = filda_dict[band_lst_emit_bg].values
	
	data_bg = BackgroundData(sig_bg, 
	                         sigma_bg, 
	                         filda_dict['Sensor_Zenith'].values, 
	                         filda_dict['FP_Area'].values,
	                         emit_bg,
	                         filda_dict['LST_Night_1KM'].values,
	                         filda_dict['LST_Night_1KM_STD'].values,
	                         filda_dict['LST_Day_1KM'].values,
	                         filda_dict['LST_Day_1KM_STD'].values)
	
	# process the fire signal
	snr_fire = np.array([nl.SNR_fire[band] for band in nl.sel_fire_bands])
	band_lst_fire = ['FP_' + band + '_Rad' for band in nl.sel_fire_bands]
	sig_fire = filda_dict[band_lst_fire].values
	
	band_lst_mean = ['FP_' + band + '_Rad_Mean' for band in nl.sel_fire_bands]
	sig_mean = filda_dict[band_lst_mean].values
	sigma_fire = (sig_fire - sig_mean) / np.reshape(snr_fire, (1,) + snr_fire.shape)
	sigma_fire = sig_fire / np.reshape(snr_fire, (1,) + snr_fire.shape)
	
	data_fire = FireData(sig_fire - sig_mean, 
	                     sigma_fire, 
	                     filda_dict['Sensor_Zenith'].values, 
	                     filda_dict['FP_Area'].values,
	                     filda_dict['FP_Power'].values,
	                     filda_dict['FP_Gas_Flaring'].values.astype(int),
	                     filda_dict['Static_flag'].values.astype(int),
	                     filda_dict['FP_Bowtie'].values,
	                     sig_fire)
	
	return data_fire, data_bg


class TensorSettings:
    def __init__(self, band_list):
        # Initialize shared tensors as attributes of the class
        self.obs   = theano.shared(20000*np.ones(len(band_list)))
        self.obs_sigma = theano.shared(np.ones(len(band_list)))
        self.vza    = theano.shared(0.)
        self.emit   = theano.shared(np.ones(len(band_list)))
        self.C      = theano.shared(0.5)
        self.C_sigma  = theano.shared(0.5)
        self.tb     = theano.shared(300.0)
        self.tb_sigma = theano.shared(5.0)
        self.frp    = theano.shared(0.2)
        self.frp_sigma = theano.shared(0.2)
        self.area = theano.shared(600000)
        self.frac = theano.shared(0.1)
        self.frac_sima = theano.shared(0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def init_background_estimation(nl, bgs):
	'''
	Initialize the background temperature estimator
	
	Parameters:
	nl: namelist object, configuration for the estimator
	bgs: background sensor object, sensor configuration
	
	'''
	
	# Start to define the background estimation model
	
    # Define some pseudo values as placeholders to define the MCMC model
    # Placeholders means it will only be used to create the computation 
    # graph, Those values will be updated for each real case

	ts_bg = TensorSettings(nl.sel_bg_bands)
	
	ts_bg.C.set_value(nl.mean_C)
	ts_bg.C_sigma.set_value(nl.sigma_C)
	
	ts_bg.tb.set_value(nl.mean_tb)
	ts_bg.tb_sigma.set_value(nl.sigma_tb)


	with pm.Model() as bg_estimator:
		# Prior for 't_b'
		t_b = pm.Normal('t_b', mu=ts_bg.tb, sigma=ts_bg.tb_sigma)
		# Prior for 'C'
		C   = pm.Normal('C',   mu=ts_bg.C,  sigma=ts_bg.C_sigma)
		
		# well, I have to use a weird line breaking, to make the line 
		# shorter for reading
		likelihood = pm.Normal(
			'tb_obs',
			# some of the tensor are used in define the dynamic graph...
			mu=MT.get_band_radiance_BG_tensor(t_b, C, ts_bg, bgs),
			sigma=ts_bg.obs_sigma,
			observed=ts_bg.obs
		)
		                       
	return bg_estimator, ts_bg

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def init_biphasic_estimation(nl, fss):
	'''
	Initialize the bi_phasic fire estimator
	
	Parameters:
	nl: namelist object, configuration for the estimator
	fss: background sensor object, sensor configuration
	
	'''
	
	# Start to define the bi_phasic estimation model
	
	# Define some pseudo values as placeholders to define the MCMC model
	# Placeholders means it will only be used to create the computation 
	# graph, Those values will be updated for each real case
	ts_bi_fire = TensorSettings(nl.sel_fire_bands)

	with pm.Model() as biphase_estimator:
		# For PyMC3
		if is_pymc3:
			if nl.flag_dist == 'G':
			
				PositiveNormal_P = pm.Bound(pm.Normal, lower=nl.p_min)

# 				PositiveNormal_Ts = pm.Bound(pm.Normal, 
# 				                             lower=300, 
# 				                             upper=2200)
# 				                             
# 				PositiveNormal_Tf = pm.Bound(pm.Normal,
# 				                             lower=300, 
# 				                             upper=2200)				
				
				PositiveNormal_Ts = pm.Bound(pm.Normal, 
				                             lower=nl.mean_ts-nl.sigma_ts, 
				                             upper=nl.mean_ts+nl.sigma_ts)
				                             
				PositiveNormal_Tf = pm.Bound(pm.Normal, 
				                             lower=nl.mean_tf-nl.sigma_tf, 
				                             upper=nl.mean_tf+nl.sigma_tf)				

				# Prior for 't_s'
				t_s = PositiveNormal_Ts('t_s', mu=nl.mean_ts, sigma=nl.sigma_ts)
				# Prior for 't_f'
				t_f = PositiveNormal_Tf('t_f', mu=nl.mean_tf, sigma=nl.sigma_tf)
				# Prior for 'p_s'  
				p_s = PositiveNormal_P('p_s', mu=nl.mean_ps, sigma=nl.sigma_ps)
				# Prior for 'p_f'
				p_f = PositiveNormal_P('p_f', mu=nl.mean_pf, sigma=nl.sigma_pf)
				
			if nl.flag_dist == 'U':
				# Prior for 't_s'
				t_s = pm.Uniform('t_s', lower=nl.mean_ts-nl.sigma_ts, 
				                        upper=nl.mean_ts+nl.sigma_ts)
				# Prior for 't_f'
				t_f = pm.Uniform('t_f', lower=nl.mean_tf-nl.sigma_tf, 
				                        upper=nl.mean_tf+nl.sigma_tf)
				# Prior for 'p_s'  
				p_s = pm.Uniform('p_s', lower=nl.p_min, 
				                        upper=nl.mean_ps+nl.sigma_ps)
				# Prior for 'p_f'
				p_f = pm.Uniform('p_f', lower=nl.p_min, 
				                        upper=nl.mean_pf+nl.sigma_pf)	
			
		else:
			if nl.flag_dist == 'G':
				t_s = pm.TruncatedNormal('t_s', mu=nl.mean_ts, sigma=nl.sigma_ts, lower=300, upper = 2200)
				# Prior for 't_f'
				t_f = pm.TruncatedNormal('t_f', mu=nl.mean_tf, sigma=nl.sigma_tf, lower=300, upper = 2200)
				# Prior for 'p_s'  
				p_s = pm.TruncatedNormal('p_s', mu=nl.mean_ps, sigma=nl.sigma_ps, lower=nl.p_min)
				# Prior for 'p_f'
				p_f = pm.TruncatedNormal('p_f', mu=nl.mean_pf, sigma=nl.sigma_pf, lower=nl.p_min)  			
		
			if nl.flag_dist == 'U':
				# Prior for 't_s'
				t_s = pm.Uniform('t_s', lower=nl.mean_ts-nl.sigma_ts, 
				                        upper=nl.mean_ts+nl.sigma_ts)
				# Prior for 't_f'
				t_f = pm.Uniform('t_f', lower=nl.mean_tf-nl.sigma_tf, 
				                        upper=nl.mean_tf+nl.sigma_tf)
				# Prior for 'p_s'  
				p_s = pm.Uniform('p_s', lower=nl.p_min, 
				                        upper=nl.mean_ps+nl.sigma_ps)
				# Prior for 'p_f'
				p_f = pm.Uniform('p_f', lower=nl.p_min, 
				                        upper=nl.mean_pf+nl.sigma_pf)
					
		# well, I have to use a weird line breaking, to make the line 
		# shorter for reading
		
		# Define the radiance likelihood
		likelihood_bi_rad = pm.Normal(
			'bi_obs', 
			mu=MT.get_band_biphasic_radiance_tensor(t_s, t_f, p_s, p_f, 
													ts_bi_fire, fss),
			sigma=ts_bi_fire.obs_sigma,
			observed=ts_bi_fire.obs
		)
		
		likelihood_bi_frp = pm.Normal(
			'bi_frp', 
			mu=MT.get_bi_frp_tensor(t_s, t_f, p_s, p_f, ts_bi_fire), 
			sigma=ts_bi_fire.frp_sigma,
			observed=ts_bi_fire.frp)
						   
	return biphase_estimator, ts_bi_fire

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def init_uniphasic_estimation(nl, fss):
	'''
	Initialize the bi_phasic fire estimator
	
	Parameters:
	nl: namelist object, configuration for the estimator
	fss: background sensor object, sensor configuration
	
	'''
	
	# Start to define the bi_phasic estimation model
	
	# Define some pseudo values as placeholders to define the MCMC model
	# Placeholders means it will only be used to create the computation 
	# graph, Those values will be updated for each real case
	printf(f'Define tensor for Uniphasic estimator...', 1, prefix = 'MCBEF')
	ts_uni_fire = TensorSettings(nl.sel_fire_bands)
	
	with pm.Model() as uniphase_estimator:

		if is_pymc3:
			if nl.flag_dist == 'G':
				PositiveNormal_P = pm.Bound(pm.Normal, lower=nl.p_min)
				
				PositiveNormal_T = pm.Bound(pm.Normal, lower=nl.mean_t-nl.sigma_t, \
				                            upper = nl.mean_t+nl.sigma_t)
				# Prior for 't_mean'
				t_mean = PositiveNormal_T('t_mean', mu=nl.mean_t, sigma=nl.sigma_t)
				# Prior for 'p_f'
				p_mean = PositiveNormal_P('p_mean', mu=nl.mean_p, sigma=nl.sigma_p)

			if nl.flag_dist == 'U':
				# Prior for 't_mean'
				t_mean = pm.Uniform('t_mean', lower=nl.mean_t-nl.sigma_t, 
				                              upper=nl.mean_t+nl.sigma_t)
				# Prior for 'p_mean'
				p_mean = pm.Uniform('p_mean', lower=nl.p_min, 
				                              upper=nl.mean_p+nl.sigma_p)

			
		else:
			if nl.flag_dist == 'G':
				t_mean = pm.TruncatedNormal('t_mean', mu=nl.mean_t, sigma=nl.sigma_t, lower=300, upper = 2200)
				# Prior for 'p_f'
				p_mean = pm.TruncatedNormal('p_mean', mu=nl.mean_p, sigma=nl.sigma_p, lower=nl.p_min)
				
			if nl.flag_dist == 'U':
				# Prior for 't_s'
				t_mean = pm.Uniform('t_mean', lower=nl.mean_t-nl.sigma_t, 
				                              upper=nl.mean_t+nl.sigma_t)
				# Prior for 't_f'
				p_mean = pm.Uniform('p_mean', lower=nl.p_min, 
				                              upper=nl.mean_p+nl.sigma_p)

		# well, I have to use a weird line breaking, to make the line 
		# shorter for reading

		# Define the radiance likelihood
		likelihood_bi_rad = pm.Normal(
			'uni_obs', 
			mu=MT.get_band_uniphasic_radiance_tensor(t_mean, p_mean, ts_uni_fire, fss),
			sigma=ts_uni_fire.obs_sigma,
			observed=ts_uni_fire.obs
		)
		
		likelihood_bi_frp = pm.Normal(
			'uni_frp', 
			mu=MT.get_uni_frp_tensor(t_mean, p_mean, ts_uni_fire), 
			sigma=ts_uni_fire.frp_sigma,
			observed=ts_uni_fire.frp)
						   
	return uniphase_estimator, ts_uni_fire


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def estimate_bg(bg_estimator, nl, ts_bg, bg_obs, verbose = 0):

	ts_bg.obs.set_value(bg_obs['obs'])
	ts_bg.obs_sigma.set_value(bg_obs['obs_sigma'])
	ts_bg.vza.set_value(bg_obs['vza'])
	ts_bg.emit.set_value(bg_obs['emit'])
	
	ts_bg.emit.set_value(bg_obs['emit'])
	ts_bg.emit.set_value(bg_obs['emit'])
	
	ts_bg.tb.set_value(bg_obs['lst'])
	ts_bg.tb_sigma.set_value(bg_obs['lst_sigma'])

	with bg_estimator:
		map_estimate = pm.find_MAP(method='L-BFGS-B', #SLSQP #'L-BFGS-B'
		                           progressbar=nl.flag_verbose)
		                           
	est_tb = map_estimate['t_b']
	est_C  = map_estimate['C']
	
	message = f"Tb: {est_tb:4.2f} K; {est_tb-273.15:4.2f} C; Scale: {est_C:4.3f}"
	printf(message, verbose, prefix = 'MCBEF')
	
	return map_estimate
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def swap_variables(data_dict, key1, key2):
	"""
	Swap the values of key1 and key2 in data_dict.
	
	Parameters:
	data_dict (dict): The dictionary containing the variables.
	key1 (str): The first key to swap.
	key2 (str): The second key to swap.
	
	Returns:
	dict: The updated data_dict with key1 and key2 values swapped.
	"""
	data_dict[key1], data_dict[key2] = data_dict[key2], data_dict[key1]
	return data_dict
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def estimate_fire(estimator, nl, ts_fire, fire_obs, verbose = 0):
	'''	
	For given fire observation and designated estimator, estimate_fire 
	updates the fire tensor and conducts samplings.
	'''

	ts_fire.obs.set_value(fire_obs['obs'])
	ts_fire.obs_sigma.set_value(fire_obs['obs_sigma'])
	ts_fire.vza.set_value(fire_obs['vza'])
	ts_fire.C.set_value(fire_obs['C'])
	ts_fire.frp.set_value(fire_obs['frp'])
	ts_fire.frp_sigma.set_value(fire_obs['frp_sigma'])
	
	if 't_s' in estimator.named_vars.keys():
		model_type = 2
		# get the variables
		t_s = estimator['t_s']
		t_f = estimator['t_f']
		p_s = estimator['p_s']
		p_f = estimator['p_f']

	else:
		model_type = 1
		
	
	# the estimator can be biphasic or uniphasic
	traces = []
	
	if nl.flag_grad_sampling:
	
		message = f"Choosing No U-Turn Sampler (NUTS)..."
		printf(message, 1, prefix = 'MCBEF')		
		with estimator:
			for i in range(nl.num_chain):
				t_ = pm.sample(draws = nl.num_draw, tune=nl.num_tune, 
							   step = pm.NUTS(target_accept=0.9),
							   return_inferencedata=False, 
							   cores = 1, 
							   chains = 1,
							   progressbar=nl.flag_verbose)
				traces.append(t_)
	else:	
		with estimator:
			map_estimate = pm.find_MAP(method='L-BFGS-B', progressbar=False)
		
			if model_type == 2 :

				if (map_estimate['t_s'] > map_estimate['t_f']):
				
					# Swap t_s and t_f
					map_estimate = swap_variables(map_estimate, 
					                              't_s', 
					                              't_f')		

					map_estimate = swap_variables(map_estimate, 
					                              'p_s', 
					                              'p_f')
					if nl.flag_dist == 'G':
						map_estimate = swap_variables(map_estimate, 
													  't_s_interval__', 
													  't_f_interval__')					
						map_estimate = swap_variables(map_estimate, 
													  'p_s_lowerbound__', 
													  'p_f_lowerbound__')
					
		for i in range(nl.num_chain):
		
			with estimator:
				
				if model_type == 2 :
					# Define custom step methods
					step_ts_ps = pm.Metropolis(vars=[t_s, p_s])
					step_tf_pf = pm.Metropolis(vars=[t_f, p_f])
					
					step = [step_ts_ps, step_tf_pf]				
				else:
					step = pm.Metropolis()
				
			
				t_ = pm.sample(draws = nl.num_draw, tune=nl.num_tune, 
							   step = step, 
							   return_inferencedata=False, 
							   cores = 1, 
							   chains = 1, 
							   start=map_estimate, 
							   progressbar=nl.flag_verbose)
				traces.append(t_)
							  
	varnames = [varname for varname in traces[0].varnames if not 
	            (varname.endswith('_lowerbound__') or '_interval__' in varname)]
	trace = {}
	for varname in varnames:
		trace[varname] = []
	
	for item in traces:
		for key in trace.keys():
			trace[key] = trace[key] + item[key].tolist()
	
	for key in trace.keys():
		trace[key] = np.array(trace[key])
# 	print(trace.keys())
	return trace

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def is_valid_data(data):
	""" 
	Check if the provided data dictionary contains NaN values. 
	"""
	for key, value in data.items():
		# Check for NaN values in the numpy array or scalar
		if np.isnan(value).any():  
			return False
	return True

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def estimate_one(nl, i, data_bg, data_fire,
                 bg_estimator, biphase_estimator, uniphase_estimator,
                 ts_bg, ts_bi_fire, ts_uni_fire, 
                 bgs, fss, verbose = 0):
	
	add_flag = 0
	
	bg_obs = {
		'obs': data_bg[i][0],
		'obs_sigma': data_bg[i][1],
		'vza': data_bg[i][2],
		'emit': data_bg[i][4],
		'lst': data_bg[i][5],
		'lst_sigma': data_bg[i][6]
	}
	
	fire_obs = {
		'obs': data_fire[i][0],
		'obs_sigma': data_fire[i][1],
		'vza': data_fire[i][2],
		'area': data_fire[i][3]
	}
	
	flag_gas_flaring = data_fire[i][5]
	flag_static = data_fire[i][6]
	bowtie = data_fire[i][7]
	
	# - - - - - - - - - - - -
	# Validate data before processing
	if bowtie>0.1:
		message = f"Bowtie detection at {i}"
		printf(message, 1, prefix = 'MCBEF')
		return FILLVALUE, FILLVALUE, FILLVALUE_FRP, FLAG_BOWTIE
	
	if not is_valid_data(bg_obs):
		message = f"Invalid background data encountered at {i}"
		printf(message, 1, prefix = 'MCBEF')
		return FILLVALUE, FILLVALUE, FILLVALUE_FRP, FLAG_MISSING_BG
	
	if not is_valid_data(fire_obs):
		message = f"Invalid fire data encountered at {i}"
		printf(message, 1, prefix = 'MCBEF')	
		return FILLVALUE, FILLVALUE, FILLVALUE_FRP, FLAG_MISSING_FIRE

	# - - - - - - - - - - - -
	# estimate the background temperature and scaling factor
	try:
		est_bg = estimate_bg(bg_estimator, nl, ts_bg, bg_obs, 
							 verbose = verbose)
	except (pm.exceptions.SamplingError, RuntimeError) as e:
		message = f" WARINING!! Background estimation failed, use climatology..."
		est_bg = {}
		est_bg['C'] = nl.mean_C
		est_bg['t_b'] = bg_obs['lst']
		printf(message, 1, prefix = 'MCBEF')
		add_flag = FLAG_BACKGROUND
		
	# update the estimated scaling factor
	fire_obs['C'] = est_bg['C']
	
	# - - - - - - - - - - - -
	# FRP correction
	frp_idx = nl.sel_fire_bands.index('M13')
	product = np.exp(-(fire_obs['C']*fss.tau_wvp+fss.tau_other_gas) \
					 / np.cos(np.deg2rad(fire_obs['vza'])) ) * fss.rsr
	tt = np.trapz(product, fss.lambdas, axis = 1)/\
		 np.trapz(fss.rsr, fss.lambdas, axis = 1)
	wl_fire = np.nanmean(fss.lambdas, axis = 1)
# 	print('wl_fire', wl_fire, tt)
	correct_rad = fire_obs['obs']/tt
	
# 	print(f'correct_rad, {correct_rad}')
	
	f_function = interpolate.interp1d(wl_fire, correct_rad)
	inter_wl = np.linspace(wl_fire[0], wl_fire[-1], 200)
	inter_rad = f_function(inter_wl) 
	irrad = np.trapz(inter_rad, inter_wl)
	# 	print('FRP area under the curve', )
	FRP_correct = irrad*data_fire[i][3]*2*np.pi*1e-6
	
	FPR_R_AC = data_fire[i][4]/tt[frp_idx]
	# 	print('FRP_correct', FRP_correct)
	
	# 	fire_obs['obs_sigma'] = fire_obs['obs_sigma']/tt
	fire_obs['frp']       = FPR_R_AC * 1.2
	# 	fire_obs['frp']       = FPR_R_AC*0.8 + FRP_correct * 0.2 # * 1.1
	fire_obs['frp_sigma'] = fire_obs['frp'] * nl.frp_sigma_scale
	

    # - - - - - - - - - - - -
	# draw samples...
	if (fire_obs['frp'] > nl.thd_frp) & (flag_gas_flaring <=0) : # & (flag_static <=0)
		
		try:
			flag_mode = FLAG_BIPHASIC
			trace = estimate_fire(biphase_estimator, nl, ts_bi_fire, fire_obs, 
		                      	  verbose = verbose)
		                      	  
		except (pm.exceptions.SamplingError, RuntimeError) as e:
			message = f"WARINING!! Sampling on biphasic attempt failed, adjusting to uniphasic..."
			printf(message, 1, prefix = 'MCBEF')
			
			try:
				flag_mode = FLAG_DEDGRADATE + add_flag
				trace = estimate_fire(uniphase_estimator, nl, ts_uni_fire, fire_obs, 
								      verbose = verbose)
								      
			except (pm.exceptions.SamplingError, RuntimeError) as e:
				flag_mode = FLAG_FAIL_BIPHASIC
				message = f"WARINING!! Fail biphasic at {i}"
				printf(message, 1, prefix = 'MCBEF')
				return FILLVALUE, FILLVALUE, FILLVALUE_FRP, flag_mode

	else:
		try:
			flag_mode = FLAG_UNIPHASIC + add_flag
			trace = estimate_fire(uniphase_estimator, nl, ts_uni_fire, fire_obs, 
							  verbose = verbose)
							  
		except (pm.exceptions.SamplingError, RuntimeError) as e:
			flag_mode = FLAG_FAIL_UNIPHASIC
			message = f"WARINING!! Fail uniphasic at {i}"
			printf(message, 1, prefix = 'MCBEF')
			return FILLVALUE, FILLVALUE, FILLVALUE_FRP, flag_mode
		
# 	est_bg = None
# 	trace = None
# 	flag_mode = None
	return est_bg, trace, FPR_R_AC, flag_mode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def estimate_batch(nl, data_bg, data_fire, 
                   bg_estimator, biphase_estimator, uniphase_estimator,
                   ts_bg, ts_bi_fire, ts_uni_fire,
                   bgs, fss, interest_pos):

	rows = interest_pos[1] - interest_pos[0]
	result = dict(x=[], est_bg = [], trace = [], time = [], status=[], FPR_R_AC=[])    
	
	for this_point in range(interest_pos[0], interest_pos[1]):
	
		start_time = time.time()
		# MZ add atmospheric corrected (AC) FRP output
		est_bg, trace, FPR_R_AC, flag_mode = estimate_one(nl, this_point, data_bg, 
												 data_fire, bg_estimator, 
												 biphase_estimator, 
												 uniphase_estimator,
												 ts_bg, ts_bi_fire, ts_uni_fire,
												 bgs, fss, 
												 verbose = nl.flag_verbose)
		end_time = time.time()  # End time measurement
		execution_time = end_time - start_time  # Calculate execution time
		print(f" - Execution time: {this_point} {execution_time:4.2f} seconds\n")
		
		result['x'].append(this_point)
		result['time'].append(execution_time)
		result['est_bg'].append(est_bg)
		result['trace'].append(trace)
		result['status'].append(flag_mode)
		result['FPR_R_AC'].append(FPR_R_AC)
	
	return result

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def MCBEF_MP(nl, data_bg, data_fire, 
             bg_estimator, biphase_estimator, uniphase_estimator,
             ts_bg, ts_bi_fire, ts_uni_fire,
             bgs, fss):

	len_array = np.arange(0, len(data_fire))
	sub_arrays_list = split_array(len_array, nl.num_core)
	
	# MZ Oct 7 2024, provide an option to automatically determine 
	# the number of core used in MP process
	if nl.num_core == -1:
		nl.num_core = (multiprocessing.cpu_count() - 4)
	
	pool = Pool(processes=nl.num_core)
	
	func = partial(estimate_batch,
				   nl, data_bg, data_fire, 
				   bg_estimator, biphase_estimator, uniphase_estimator,
				   ts_bg, ts_bi_fire, ts_uni_fire,
				   bgs, fss)            
			 
	result = pool.map(func, sub_arrays_list)
	
	pool.close()
	pool.join()
	
	output = {}
	for item in result:
		for i, tr, ti, bg, frp, st in zip(item['x'],    item['trace'], 
									      item['time'], item['est_bg'], 
									      item['FPR_R_AC'], item['status']):
			output[i] = {}
			output[i]['trace'] = tr
			output[i]['time'] = ti
			output[i]['est_bg'] = bg
			output[i]['FP_Power_R_AC'] = frp
			output[i]['status'] = st
	
	return output

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_sample_output(array_dict, filda_dict):
	'''
	'''
	# prepare the sample output
	ops = set(filda_dict['overpass'])
	
	sample_dicts = {}
	for op in ops:
		sample_dicts[op] = {}
		idx = np.where(filda_dict['overpass']== op)
		for key in array_dict.keys():
			if array_dict[key].ndim == 1:
				sample_dicts[op][key] = array_dict[key][idx[0]]
			if array_dict[key].ndim == 2:
				sample_dicts[op][key] = array_dict[key][idx[0], :]
	
		params = ['FP_confidence', 'FP_Land_Type', 'FP_Gas_Flaring', 'Static_flag',
				  'FP_Peatland', 'FP_Peatfrac', 'FP_SAA_flag', 'FP_Latitude',
				  'FP_Longitude','FP_Area', 'FP_Line', 'FP_Sample', 'FP_Bowtie']
		for param in params:
			sample_dicts[op][param] = filda_dict[param].values[idx[0]]
		sample_dicts[op]['FP_Power_R'] = filda_dict['FP_Power'].values[idx[0]]
	
	return sample_dicts
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def get_state_output(array_dict, filda_dict):
	'''
	'''
	n_fire = array_dict['t_s'].shape[0]
	
	state_dict = {}
	state_params = ['t_s', 'p_s', 't_f', 'p_f', 'FP_Power_T', 'FP_Power_F', 'FP_Power_S']
	for param in state_params:
		state_dict[param + '_low']  = np.full(n_fire, np.nan)
		state_dict[param + '_upp']  = np.full(n_fire, np.nan)
		state_dict[param + '_mean'] = np.full(n_fire, np.nan)
		state_dict[param + '_mode'] = np.full(n_fire, np.nan)
		state_dict[param + '_sd']   = np.full(n_fire, np.nan)
	
	# make sure QA is int...
	_QA = array_dict['QA_flag'].astype(int)
	for i in range(n_fire):
		
		QA_flag = _QA[i]
		t_b = array_dict['t_b'][i]
		# in the case that uniphasic model is used
		if (QA_flag == FLAG_UNIPHASIC) | (QA_flag == FLAG_DEDGRADATE) | \
			(QA_flag == int(FLAG_UNIPHASIC+FLAG_BACKGROUND)) | \
			(QA_flag == int(FLAG_DEDGRADATE+FLAG_BACKGROUND)):
	
			# we do put enough constrain to make sure the samplings are
			# within a physics explainable range, however, let's still
			# put this check to ensure the statics summaries are generated
			# based on a set of reasonable data 
			valid = np.where( (array_dict['t_f'][i, :]>t_b) & \
							  (array_dict['p_f'][i, :]>0) )
	
			state_params = ['t_f', 'p_f', 'FP_Power_T', 'FP_Power_F']  
			for param in state_params:
				x = array_dict[param][i]
				x = x[valid]
				hdi = pm.hdi(x, hdi_prob=0.95)
				state_dict[param + '_low'][i] = hdi[0]
				state_dict[param + '_upp'][i] = hdi[1]
	
				state_dict[param + '_mean'][i] = np.nanmean(x)
				state_dict[param + '_sd'][i]   = np.nanstd(x, ddof = 1)
				state_dict[param + '_mode'][i] = stats.mode(x)[0]
	
		if (QA_flag == FLAG_BIPHASIC) | \
			(QA_flag == int(FLAG_BIPHASIC+FLAG_BACKGROUND)) :
	
			valid = np.where( (array_dict['t_s'][i, :]>t_b) & \
							  (array_dict['p_s'][i, :]>0)   & \
							  (array_dict['t_f'][i, :]>t_b) & \
							  (array_dict['p_f'][i, :]>0))
	
			state_params = ['t_s', 'p_s', 't_f', 'p_f', 'FP_Power_T', 'FP_Power_F', 'FP_Power_S']
			for param in state_params:
				x = array_dict[param][i]
				x = x[valid]
				hdi = pm.hdi(x, hdi_prob=0.95)
				state_dict[param + '_low'][i] = hdi[0]
				state_dict[param + '_upp'][i] = hdi[1]
	
				state_dict[param + '_mean'][i] = np.nanmean(x)
				state_dict[param + '_sd'][i]   = np.nanstd(x, ddof = 1)
				state_dict[param + '_mode'][i] = stats.mode(x)[0]				
	
	for key in array_dict.keys():
		if array_dict[key].ndim == 1:
			state_dict[key] = array_dict[key]
	
	params = ['FP_confidence', 'FP_Land_Type', 'FP_Gas_Flaring', 
	          'Static_flag', 'FP_Peatland', 'FP_Peatfrac', 
	          'FP_SAA_flag', 'FP_Latitude', 'FP_Longitude','FP_Area', 
	          'FP_Line', 'FP_Sample', 'FP_Bowtie']
	for param in params:
		state_dict[param] = filda_dict[param].values
		 
	state_dict['FP_Power_R'] = filda_dict['FP_Power'].values
	state_dict['overpass'] = ST.convert_to_interval_index(filda_dict['overpass'].values)
	
	FP_combined_land_type = copy.deepcopy(state_dict['FP_Land_Type'])
	
	idx = np.where(state_dict['Static_flag'] == 1)
	FP_combined_land_type[idx] = 101
	
	idx = np.where(state_dict['FP_Gas_Flaring'] == 1)
	FP_combined_land_type[idx] = 102
	
	state_dict['FP_combined_land_type'] = FP_combined_land_type
	
	return state_dict				


def post_processing(filda_dict, result, nl):
	'''
	'''    
	n_fire = len(result)
	n_samples = nl.num_chain*nl.num_draw
	
	
	array_dict = {}
	array_dict['t_b']       = np.full(n_fire, np.nan)
	array_dict['C']         = np.full(n_fire, np.nan)
	array_dict['QA_flag']   = np.full(n_fire, np.nan)
	array_dict['FP_Power_R_AC']   = np.full(n_fire, np.nan)
	array_dict['t_s']        = np.full( (n_fire, n_samples), np.nan)
	array_dict['p_s']        = np.full( (n_fire, n_samples), np.nan)
	array_dict['t_f']        = np.full( (n_fire, n_samples), np.nan)
	array_dict['p_f']        = np.full( (n_fire, n_samples), np.nan)
	array_dict['FP_Power_T'] = np.full( (n_fire, n_samples), np.nan)
	array_dict['FP_Power_F'] = np.full( (n_fire, n_samples), np.nan)
	array_dict['FP_Power_S'] = np.full( (n_fire, n_samples), np.nan)
	
	sigma = 5.6704*10**-8
	
	for i in range(len(result)):
	
		qa = int(result[i]['status'])
		array_dict['QA_flag'][i] = qa
		array_dict['FP_Power_R_AC'][i] = result[i]['FP_Power_R_AC']
		area = filda_dict['FP_Area'].values[i]
		
		if qa >= 100:
			continue
		else:
			array_dict['t_b'][i] = result[i]['est_bg']['t_b']
			array_dict['C'][i]   = result[i]['est_bg']['C']
	
			if (qa == FLAG_UNIPHASIC) | (qa == FLAG_DEDGRADATE) | \
				(qa == int(FLAG_UNIPHASIC+FLAG_BACKGROUND))| \
				(qa == int(FLAG_DEDGRADATE+FLAG_BACKGROUND)):
				array_dict['t_f'][i, :] = result[i]['trace']['t_mean']
				array_dict['p_f'][i, :] = result[i]['trace']['p_mean']
				
				array_dict['FP_Power_T'][i, :] = area *sigma* \
				result[i]['trace']['t_mean']**4 * result[i]['trace']['p_mean'] * 1e-6

				array_dict['FP_Power_F'][i, :] = area *sigma* \
				result[i]['trace']['t_mean']**4 * result[i]['trace']['p_mean'] * 1e-6

			if (qa == FLAG_BIPHASIC) | (qa == int(FLAG_BIPHASIC+FLAG_BACKGROUND)):
				array_dict['t_s'][i, :] = result[i]['trace']['t_s']
				array_dict['p_s'][i, :] = result[i]['trace']['p_s'] 
				array_dict['t_f'][i, :] = result[i]['trace']['t_f']
				array_dict['p_f'][i, :] = result[i]['trace']['p_f']
				
				array_dict['FP_Power_T'][i, :] = area *sigma* \
				(result[i]['trace']['t_s']**4 * result[i]['trace']['p_s'] + \
				 result[i]['trace']['t_f']**4 * result[i]['trace']['p_f']) * 1e-6
				 
				array_dict['FP_Power_F'][i, :] = area *sigma* \
				(result[i]['trace']['t_f']**4 * result[i]['trace']['p_f']) * 1e-6
				
				array_dict['FP_Power_S'][i, :] = area *sigma* \
				(result[i]['trace']['t_s']**4 * result[i]['trace']['p_s']) * 1e-6				
	  
	sample_dicts = get_sample_output(array_dict, filda_dict)
				
	state_dict   = get_state_output(array_dict, filda_dict)
	
	return sample_dicts, state_dict
