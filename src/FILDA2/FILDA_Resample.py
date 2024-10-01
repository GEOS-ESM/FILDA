#---------------------------------------
# function of resampling DNB to M band
# Comment in progress
#---------------------------------------


import numpy as np
from netCDF4 import Dataset

#----------------------------------------------------------------------------
def resample_DNB2MOD(dataset, lut, invalid=None, var=['DNB_observations']):	

	'''
	resample_DNB2MOD resample the variables from DNB grid to M-band

	-----------
	Parameters:
	dataset		: directory for holding the data in DNB grids
	lut			: DNB to M band resampling lookup table
	
	-----------
	Returns:
	output		: directory for holding resample results in M-band grid
	
	'''
	
	print(' - Resampling DNB data to M band...\n')
	
	output = {}
	for key in var:
		output[key], overlapIdx = resample_DNB_to_Mband(dataset[key], lut)
	
	output['Overlapped_Indices'] = overlapIdx
	# make the negative value into NaN
	output['DNB_observations'][np.where(output['DNB_observations']<=0)] = np.nan
	# if M band invalid pixel index is provided, use it to further mask out
	# invalid pixel
	if invalid != None:
		output['DNB_observations'][invalid]	= np.nan
	
	return output 

#----------------------------------------------------------------------------
def resample_DNB_to_Mband(observation, reLUT):
	'''
	resample_DNB_to_Mband resamples one variables from DNB to M band

	-----------
	Parameters:
	observation : array of the data in DNB grids
	lut			: DNB to M band resampling lookup table
	
	-----------
	Returns:
	Granule_Weighted_Values	: Array of the resampled dataset
	Overlapped_Indices		: Array for recording the overlap ratio for each
							  pixel
	
	'''
	
	lut = Dataset(reLUT, 'r')
	Matched_DNBIDX = lut.variables['Matched_DNBIDX'][:]
	Weights        = lut.variables['Weights'][:]	

	DNB_Radiance = observation

	# the original weight has sum(weights) > 0, 1 > sum(weights) > 0
	# here we modify the weights to keep the sum of the weight to 1
	# print the sum_W to see different conditions
# 	sum_W = np.nansum(Weights, axis = 0)
# 	idx = np.where( (sum_W >0) )
# 	Weights[:, idx] = Weights[:, idx] / sum_W[idx]
	
	# Since DNB may have some negative values which has no physical meaning
	# we need to use DNB to reconstruct the weight again...
	
	# Number of columns for MBand granule
	Granule_Columns = 3200
	
	# Number of rows for MBand granule, DNB and M band have the same
	# number of rows...
	Granule_Rows = np.shape(DNB_Radiance)[0]

	# Number of rows for matched DNB pixels. The maximum nuber of DNB 
	# pixels that matched to an MBand pixel is 12. As a result, 
	# the Matched_Rows will be 12
	Matched_Rows = np.shape(Matched_DNBIDX)[0]

	#The number of elements (16*3200)
	#of a 1d MBnad array for a scan.
	Matched_Cols= np.shape(Matched_DNBIDX)[1]	

	#An aray for sorting weighted values for each scan
	Scan_Weighted_Values = np.zeros((Matched_Rows,Matched_Cols))	

	#An aray for sorting weighted values for whole granule
	Granule_Weighted_Values = np.zeros((Granule_Rows,Granule_Columns))

	#Overlapped_Indices
	Overlapped_Indices = lut.variables['MBand_Overlapped'][:]
	Overlapped_Indices = Overlapped_Indices[0:Granule_Rows,:]
	
	for k in range(0, Granule_Rows,16):
		
# 		print(' - Resampling ' + str(k//16) + ' scan...')
		
		Radiance = np.ravel(DNB_Radiance[k:k+16,:])	

		ones	      = np.ones_like(Radiance)
		idx 		  = np.where(Radiance<0)
		ones[idx] 	  = np.nan
		Radiance[idx] = np.nan
		#Finding DNB indices and weights for each MBand pixel in a scan
		Weights_mod   = np.zeros_like(Weights)
		# 0:12
		for i in range(0,Matched_Rows):
			# pass the weighted value to the variable holder
			Scan_Weighted_Values[i,:] = Radiance[Matched_DNBIDX[i,:]] * Weights[i,:]
			# also pass the modified weights to the variable holder
			Weights_mod[i,:]		  = ones[Matched_DNBIDX[i,:]] * Weights[i,:]

		# Calculating the weighted value		
		Scan_Weighted_Values_New = np.nansum(Scan_Weighted_Values,axis=0)			
		# Calculating the modified weights
		Weights_mod = np.nansum(Weights_mod, axis = 0)
		
		Scan_Weighted_Values_New = Scan_Weighted_Values_New / Weights_mod
		
		Scan_Weighted_Values_New = np.reshape(Scan_Weighted_Values_New,(16,3200))
		
		Granule_Weighted_Values[k:k+16,:] = Scan_Weighted_Values_New[:,:]
		
	Granule_Weighted_Values[Granule_Weighted_Values!=Granule_Weighted_Values] = -999.9
	
	lut.close()
	
	return Granule_Weighted_Values, Overlapped_Indices
