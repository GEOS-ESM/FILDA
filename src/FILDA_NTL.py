import glob
import numpy as np
import os

import FILDA_IO			# module for reading and writing files
import FILDA_Time_Cord  # module for manipulating the time and coordinates
import FILDA_Resample 	# modele for physical resampling from DNB to M/I band
import FILDA_Cloud 		# modele for cloud masking
import FILDA_CLT	    # modele for generating the surface light climatology


'''
data_Dir = '/Dedicated/jwang-data/shared_satData/VIIRS/yearly_data/VNP46A1/'
output_dir = './'

year_month = '201908'
tile = 'h10v05'
gen_NTL_climatology( year_month, tile, data_Dir, output_dir)

'''

def update_params_for_collection(filename, params):
	# MZ, Sept. 04
	# Make change to accommodate the changes in the C2 VNP46A1

	basename = os.path.basename(filename)
	collection = basename.split('.')[3]

	if collection == '001':
		datafiled = 'HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/'
		if "DNB_At_Sensor_Radiance" in params:
			idx = params.index("DNB_At_Sensor_Radiance")
			params[idx] = "DNB_At_Sensor_Radiance_500m"
	else:
		datafiled = 'HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields/'
		if "DNB_At_Sensor_Radiance_500m" in params:
			idx = params.index("DNB_At_Sensor_Radiance_500m")
			params[idx] = "DNB_At_Sensor_Radiance"
	
	return datafiled, params, collection
    
#-----------------------------------------------------------------------
def read_VNP46A1(filename, params):
	import h5py
	import numpy as np
	output = {}

	datafiled, params, collection = update_params_for_collection(filename, params)

	print(f' - Reading {filename}, Collection {collection}...')
	f = h5py.File(filename, 'r')

	for param in params:

		if 'QF' not in param:
		
			obs_value = f[datafiled + param][:] * 1.0

			# Check for the presence of attributes before accessing
			scale_factor = f[datafiled + param].attrs.get('scale_factor', [1.0])[0]
			add_offset = f[datafiled + param].attrs.get('add_offset', [0.0])[0]
			valid_max = f[datafiled + param].attrs.get('valid_max', [np.inf])
			valid_min = f[datafiled + param].attrs.get('valid_min', [-np.inf])
			units = f[datafiled + param].attrs.get('units', '')
			long_name = f[datafiled + param].attrs.get('long_name', '').decode('UTF-8')

			# Handle fill values
			fill_value = f[datafiled + param].attrs.get('_FillValue', None)
			if fill_value is not None:
				obs_value[obs_value == fill_value] = np.nan

			obs_value[np.where((obs_value > valid_max) | (obs_value < valid_min))] = np.nan
			obs_value = obs_value * scale_factor + add_offset
			output[param] = obs_value
			output[param + '_units'] = units
			output[param + '_long_name'] = long_name

		if 'QF' in param:
			output[param] = f[datafiled + param][:]

	# Check if "DNB_At_Sensor_Radiance" is a key in the dictionary, and replace it with "DNB_At_Sensor_Radiance_500m"
	if "DNB_At_Sensor_Radiance" in output:
		output["DNB_At_Sensor_Radiance_500m"] = output.pop("DNB_At_Sensor_Radiance")
		output["DNB_At_Sensor_Radiance_500m_units"] = output.pop("DNB_At_Sensor_Radiance_units")
		output["DNB_At_Sensor_Radiance_500m_long_name"] = output.pop("DNB_At_Sensor_Radiance_long_name")

	f.close()

	return output



#-----------------------------------------------------------------------
def cloudMask_decode(cloudMask):
	import numpy as np
	# Day night status 0=Night, 1=Day 
	day_nig_status = cloudMask & 1
	
	# Land/Water 0=Land & Desert, 1=Land & no Desert, 2=Inland Water, 3=Sea Water, 5=Coastal
	land_type = cloudMask & 14
	land_type = land_type >>1
	
	# Cloud Mask Quality 0=Poor, 1=Low, 2=Medium, 3=High
	cloudMaskQA = cloudMask & 48
	cloudMaskQA = cloudMaskQA>>4
	
	# Cloud Confidence 0=Confident Clear, 1=Probably Clear, 2=Probably Cloudy, 3=Confident Cloudy
	cloudMaskCF = cloudMask & 192
	cloudMaskCF = cloudMaskCF >>6
	
	finalCM = np.full(np.shape(cloudMask), np.nan)
	
	finalCM[np.where( (day_nig_status ==0) & (cloudMaskCF < 2) & (cloudMaskQA > 0))] = 1
	
	return finalCM


def special_handling(year_month, tile, sat):

	# MZ Sept. 05 2023, Force to use the May 2012, climatology to represent that of January to April 2012 for VNP
	# Force to use the May 2018, climatology to represent that of January to April 2012 for VJ1
	
	# Dictionary to handle special satellite and year_month conditions
	satellite_conditions = {
		'VNP': ['201201', '201202', '201203', '201204', '201205', 2012],
		'VJ1': ['201801', '201802', '201803', '201804', '201805', 2018]
	}

	if sat in satellite_conditions:
		conditions = satellite_conditions[sat]
		if year_month in conditions[:-2]:
			year_month = conditions[-2]
			print(f' - Request to create climatology for initial month of {sat}, force the year and month to {year_month}')

	# MZ Sept. 12 2023, add treatment to due the long day effects over the high latitude region
	# h04v02    
	vid = int(tile[4:])
	
	if vid < 3:
		year = int(year_month[:4])
		if sat in ['VNP', 'VJ1']:
			cutoff_year = satellite_conditions[sat][-1]
			year = cutoff_year if year <= cutoff_year else year - 1

			month = 12
			year_month = str(year) + str(month)
			print(f' - FILDA2: found High latitude region in the north hemisphere, use Sept., Oct., and Nov.')

	return year_month


def gen_NTL_climatology( year_month, tile, data_Dir, output_dir, savename, **kwargs):
	'''
	MZ Sept. 05 2023, add some modification to address the data problem for the beginning of the VIIRS mission
	'''
	number_of_day = kwargs.get('number_of_day', 90)
	lag_of_day = kwargs.get('lag_of_day', 7)
	valid_fraction = kwargs.get('valid_fraction', 0.8)
	sat = kwargs.get('sat', 'VNP')
	

	print(' - FILDA2: Generating Nighttime light climatology for', year_month)

	year_month = special_handling(year_month, tile, sat)

	input_date = year_month + '01'
	
	jdn_end = FILDA_Time_Cord.JulianDay(input_date) - lag_of_day
	jdn_ini = jdn_end - number_of_day

		
	jdns = FILDA_Time_Cord.get_date_series( FILDA_Time_Cord.GregorianDay(jdn_ini, outputformat = 'yyyymmdd'), 
											FILDA_Time_Cord.GregorianDay(jdn_end, outputformat = 'yyyymmdd'), outtype = 'nasa' )


	vnp46_var = ['DNB_At_Sensor_Radiance_500m', 'QF_DNB',  'QF_VIIRS_M12', 
				 'QF_VIIRS_M13', 'QF_VIIRS_M15', 'QF_VIIRS_M16', 'QF_Cloud_Mask']
			 
	QFvar     = ['QF_DNB',  'QF_VIIRS_M12', 'QF_VIIRS_M13', 'QF_VIIRS_M15',
				 'QF_VIIRS_M16']

	saveVars  = ['DNB_At_Sensor_Radiance_500m']

	outputs = {}
	for var in saveVars:
		outputs[var] = []
	
	filenames = []
	for jdn in jdns:

		year = jdn[1:5]
		doy  = jdn[5:]
		file_id = glob.glob(data_Dir + year + '/' + doy + '/' + '*' +  tile + '*.h5')
		if len(file_id)==0:
			print(' - Cannot find files for:', jdn, 'tile:', tile )
		else:
			filenames.append( file_id[0] )

	# MZ Sept. 05 2023, Fix the logical bug for generating the NTL...
	if len(filenames) < len(jdns)*(1-valid_fraction):
		print(f' - Can not find enough data to generate the nighttime light climatology')
		return
		
	# It will excuse only when enough files are found
	for filename in filenames:
		vnp46_data = read_VNP46A1(filename, vnp46_var)
		# decode the cloud mask
		clm = cloudMask_decode(vnp46_data["QF_Cloud_Mask"]) 
	
		# get the quality flag of the observation...
		qf  = np.full(np.shape(clm), 1.0)
	
		for var in QFvar:
			temp_qf = vnp46_data[var]
			temp_qf[np.where(temp_qf>0)] = 1
			temp_qf = 1 - temp_qf
			qf = qf * temp_qf	
	
		qf[np.where(qf<1)] = np.nan
	
		for var in saveVars:
			temp_obs = vnp46_data[var]
			temp_obs = temp_obs * qf * clm
			outputs[var].append(temp_obs)


	# calculate the statics...
	save_dict = {}
	for var in saveVars:
	
		nDay, nRow, nCol = np.shape(outputs[var])

		tempValue = np.full((nDay, nRow, nCol), np.nan)
	
		for i in range(nDay):
			tempValue[i,:,:] = outputs[var][i]	
		# do a QA...
		tempValue[np.where(tempValue!=tempValue)] = -999.
		tempValue[np.where(tempValue<=0)] = np.nan
	
		# get the mean and standard deviation of the data...	
		tempMean = np.nanmean(tempValue, axis = 0)
		tempSTD = np.nanstd(tempValue, axis = 0, ddof = 1)		
	
		save_dict[var + '_mean']  = tempMean		
		save_dict[var + '_std']   = tempSTD
		save_dict[var + '_mean' + '_units'] = vnp46_data[var + '_units']
		save_dict[var + '_std' + '_units'] = vnp46_data[var + '_units']
		save_dict[var + '_std' + '_long_name'] = 'Standard deviaion of ' + vnp46_data[var + '_long_name']
		save_dict[var + '_mean' + '_long_name']= 'Mean of ' + vnp46_data[var + '_long_name']

		if var == 'DNB_At_Sensor_Radiance_500m':

			numNone = np.full_like(tempValue, np.nan)
			numNone[np.where(tempValue == tempValue)] = 1.
			denominator = np.nansum(numNone, axis = 0)
			D = np.log(tempMean) - np.nansum(np.log(tempValue), axis = 0)/denominator
		
			invalidIdx = np.where(denominator<1)
			alpha_hat = ( 1 + (1+4*D/3)**0.5 ) / 4 / D
			beta_hat  = tempMean/alpha_hat			
		
			alpha_hat[invalidIdx] = np.nan
			beta_hat[invalidIdx] = np.nan
		
			save_dict[var + '_alpha'] = alpha_hat
			save_dict[var + '_beta']  = beta_hat
			save_dict[var + '_mean' + '_units'] = vnp46_data[var + '_units']
			save_dict[var + '_std' + '_units']  = vnp46_data[var + '_units']
			save_dict[var + '_std' + '_long_name'] = 'Standard deviaion of ' + vnp46_data[var + '_long_name']
			save_dict[var + '_mean' + '_long_name']= 'Mean of ' + vnp46_data[var + '_long_name']

		save_dict['latitude'], save_dict['longitude'] = FILDA_Time_Cord.cal_PlateCarree_grid(tile, 2400)

		FILDA_IO.write_nc_NTL( output_dir + savename, save_dict)
		
	return

























