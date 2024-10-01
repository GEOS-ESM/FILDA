
def make_namelist(): 
	data 		              = {}
	data['INPUT_DIR']         = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/LEV1B/'
	data['GEOS_DIR']          = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/LEV1B/'
	data['NTL_DIR']           = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/NIGHTTIME_LIGHT/THREE_MONTH_AVG/'
	data['SURF_DIR']          = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/LAND_COVER/'
	data['GASFLARE_DIR']      = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/GASFLARING/'
	data['PEATLAND_DIR']      = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/DATA/PEATLAND/'
	data['LUT_DIR']           = '/Dedicated/jwang-data/shared_satData/OPNL_FILDA/LUT/'
	data['OUT_DIR']           = './OUTPUT/'
	data['FIG_DIR']           = './FIG/'
	# related to the selection of fire candidates
	data['thres_BTI04']       = 295.
	data['thres_BTD_IMG']     = 10.
	data['thres_DNB']         = 0.95
	data['thres_BTD_MOD']     = 0.
	data['thres_ABSI04']      = 320
	data['thres_SATI04']      = 367
	data['thres_FLDI04']      = 208
	data['thres_FLDI05']      = 310
	data['thres_SATI05']      = 335
	# related to cloud mask
	data['thres_cloud_gross'] = 9
	data['thres_cloud_high']  = 4.5
	data['thres_cloud_water'] = 2
	data['thres_cloud_I04']	  = 295
	data['thres_cloud_I05']	  = 265
	# resolution of the temperature field...
	data['resol_lat']		  = 0.25
	data['resol_lon']		  = 0.3125
	# resolution of the contexture test
	data['thres_num']		  = 10
	data['thres_frac']		  = 0.25
	data['half_win_ini'] 	  = 5
	data['half_win_max'] 	  = 31
	data['win_step']          = 1
	# resolution of the FRP calculation
	data['thres_num_FRP']	  = 10
	data['thres_frac_FRP']	  = 0.25
	data['half_win_ini_FRP']  = 3
	data['half_win_max_FRP']  = 50
	data['win_step_FRP']      = 1


	# df = pd.DataFrame.from_dict(data, 'columns')
	# df.to_csv('namelist.input', index=False, na_rep='N/A')

	with open('namelist.input', 'w') as f:
		print(' - Writing namelist.input')

	for key in data.keys():
		with open('namelist.input', 'a') as f:
			print(key.ljust(20) +  ': '+ str( data[key] ), file = f)
