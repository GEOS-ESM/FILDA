'''
Comment in progress
'''
from scipy.optimize import minimize
from scipy.optimize import Bounds
import numpy as np 

###from pylib.firelib.FILDA_BT import *
import FILDA_BT

# lamda = dict(lamda_M07 = 0.865,
# 		 	lamda_M08 = 1.24,
# 		 	lamda_M10 = 1.61,
# 		 	lamda_M11 = 2.25,
# 		 	lamda_M12 = 3.70,
# 		 	lamda_M13 = 4.05,
# 		 	lamda_M14 = 8.55,
# 			lamda_M15 = 10.76,
# 			lamda_M16 = 12.0125,
# 	 	 	lamda_I04 = 3.757557,	#3.757494
# 		 	lamda_I05 = 11.57841)	#11.587094
		 	
def opt_func(para, lamda_all, valid_trans, measured_rad):
	### import numpy as np
	fire_T = para[0]
	bg_T   = para[1]
	frac   = para[2]

	total_rad = (frac * FILDA_BT.planck_func(lamda_all, fire_T) + (1-frac) *  FILDA_BT.planck_func(lamda_all, bg_T))*1.#*valid_trans

	return np.sum((total_rad*1e6 - measured_rad*1e6)**2)



def fire_fitting(fire_pixel):

	bands     = ['M11', 'M13', 'M14', 'M15', 'M16'] #M07, M08, M10
	lamdas    = np.array([FILDA_BT.lamda['lamda_' + band] for band in bands])
	rads      = np.concatenate( [ np.expand_dims(fire_pixel['FP_' + band + '_Rad'], axis =1) for band in bands], axis = 1)
	trans	  = np.array([FILDA_BT.transmittance[band] for band in bands])
	
	bg_temp          = []
	fire_temp        = []
	fire_frac        = []
	opt_status       = []
	
	i = 0
	for rad in rads:
		i = i+1
		valid_idx = np.where((rad==rad))
		valid_rad   = rad[valid_idx]
		valid_lamda = lamdas[valid_idx]
		valid_trans = trans[valid_idx]
		
		
		g0  = [750, 288, 0.002]
		bounds = ([550, 2200], [230, 335], [0.00005, 0.4])
		res = minimize(opt_func, g0, args=(valid_lamda, valid_trans, valid_rad), 
					   method="L-BFGS-B", bounds = bounds, tol = 1e-14,
					   options={'disp': False})
			   
		# if res.success != True:
		'''
		print('\n=============================================')
		print(i)
		print(' - FILDA: Iteration status', res.success)
		print(' - FILDA: Iteration message', res.message)
		print(' - FILDA: Opt status, %i, number of iteration: %i'%(res.status, res.nit))
		print(' - FILDA: Fire temperature: %6.4f, background temperature: %6.4f, Fire fraction: %8.6f'%(res.x[0],res.x[1], res.x[2]))
		print(' - ', valid_rad)
		'''
		
		bg_temp.append(res.x[1])
		fire_temp.append(res.x[0])
		fire_frac.append(res.x[2])
		opt_status.append(res.status)
		
	fire_pixel['FP_BG_Temp']   = np.array(bg_temp)
	fire_pixel['FP_Fire_Temp'] = np.array(fire_temp)
	fire_pixel['FP_Fire_Frac'] = np.array(fire_frac)
	fire_pixel['FP_Opt_Status']   = np.array(opt_status)

	return fire_pixel
