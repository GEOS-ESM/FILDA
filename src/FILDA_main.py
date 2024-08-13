'''

More decoration in progress...

Comment in progress

'''
import os
import sys
import numpy as np


# import FILDA
import FILDA as FI
import FILDA_Time_Cord 
import FILDA_IO
import FILDA_Resample
import FILDA_BT
import FILDA_Cloud
import FILDA_fitting


#------------------------
# User specify...
#------------------------
# jdn 		= 'A2022159' #'A2019162' 
# overpass 	= '0436' 
# sat  	    = 'VNP'

jdn         = sys.argv[1]
overpass 	= sys.argv[2]
sat 		= sys.argv[3]
print('sys.argv a python script has been called', sys.argv)

#------------------------
# main code start here
#------------------------


nameListDir = './'
flag_fires = False

# make_namelist()

# initialized the detection...
# read the necessary threshold, directory...
namelist = FI.init_detection(nameListDir)
namelist['platform'] = sat

# initialized the time
# create the necessary time string needed in the detection
time =  FILDA_Time_Cord.init_time(jdn, overpass)

#-------------------------------
# get all the files names...
#-------------------------------
file_dict = FILDA_IO.get_files(sat, namelist, time, mode = int(namelist['Layout_MODE'])) 

if file_dict:
	if file_dict['DayNightFlag'] == 'Day':
		print(' - FILDA: Skip, the daytime file was found')

	else:
		#-------------------------------
		# read DNB band data...
		#-------------------------------
		dnbData = FILDA_IO.READ_DNB(file_dict[sat+'02DNB'], file_dict[sat+'03DNB'], namelist) 

		#-------------------------------
		# read MOD band data...
		#-------------------------------
		modData, invalidIdx_mod = FILDA_IO.READ_MOD(file_dict[sat+'02MOD'], file_dict[sat+'03MOD'], namelist)

		#-------------------------------
		# read I band data...
		#-------------------------------
		imgData = FILDA_IO.READ_IMG(file_dict[sat+'02IMG'], file_dict[sat+'03IMG'], namelist)
		
		if int(namelist['COPY_STATIC']) == 1:
			namelist = FILDA_IO.copy_static(imgData, namelist, time)

		#-------------------------------
		# read geos-fp surface temperature
		#-------------------------------
		tempField = FILDA_IO.read_Geos_FP(file_dict['GEOS-FP'], ['TLML'])
		# MZ, modified the spatial resolution of the GEOS data
		namelist['resol_lat'] = tempField['resol_lat']
		namelist['resol_lon'] = tempField['resol_lon']

		#-------------------------------
		# resample DNB to M band
		#-------------------------------
		# then resample DNB radiance to M band
		DNB2MODLUT = namelist['LUT_DIR'] + sat + '_DNB2MOD_Resampling_Lookup_Table.nc'
		modData['DNB_observations'] = FILDA_Resample.resample_DNB2MOD(dnbData, DNB2MODLUT, invalidIdx_mod)['DNB_observations']
		
		#---------------------------------------
		# Calculate the brightness temperature
		#---------------------------------------
		# will remove radiance here...
		modData, imgData = FILDA_BT.get_bt(modData, imgData)

		#-------------------------------
		# M band cloud mask
		#-------------------------------
		modData, imgData = FILDA_Cloud.cloud_test(modData, imgData, tempField['TLML'][0,:,:], namelist)

		#--------------------------------------------------------------------
		# Select the candidate fire pixel and generate the background dataset
		#--------------------------------------------------------------------
		cdt_fire, posDNB = FI.sel_candidates(modData, imgData, namelist, time) 


		if cdt_fire == False:
			print(' - FILDA: No fire is found')
	
		else:	
			#-------------------------------------
			# Get the background level 1B dataset 
			#-------------------------------------
			bg_BT = FI.get_BG_IMG(imgData, modData, cdt_fire)
	
			#-------------------------------------------------
			# here we begin to identify all the fire pixels...
			# this part should be paralle computing
			#-------------------------------------------------
			fire_pixel, cdt_fire = FI.fire_test(cdt_fire, bg_BT, namelist)

			if len(fire_pixel['FP_line_mod']) <=0 :
				print(' - FILDA: No fire is found')		
			else:
				#-----------------------------------------------------
				# get the background radiance for M13 to calculate FRP
				#-----------------------------------------------------
				rad_BG, fire_pixel 	 = FI.get_BG_MOD(modData, fire_pixel)
								
				#---------------------------------------------
				# get the background radiance for M13 to 
				# calculate FRP on M band for I band detection
				#---------------------------------------------
				fire_pixel = FI.get_fire_rad13(fire_pixel, rad_BG, namelist, modData)

				#-------------------------------------------------------------
				# calculate the fire parameters on M band and split to I band
				#-------------------------------------------------------------
				MODAREALUT = namelist['LUT_DIR'] + sat + '_MOD_Pixelareas_Lookup_Table.nc'
				fire_pixel = FI.get_fire_paras(fire_pixel, MODAREALUT, namelist)
		
				
				
				if len(fire_pixel['FP_Power_QA']) <=0:
					print(' - FILDA: Skip, No fire pixel is found...')
				else:
					#-------------------------------------------------------------
					# Add the out I band information into M band output...
					#-------------------------------------------------------------
		
					fire_pixel = FILDA_fitting.fire_fitting(fire_pixel)
					fire_img, fire_mod 	= FI.add_aux_infor(fire_pixel, cdt_fire, namelist, time, modData, imgData) 

					#-------------------------------------------------------------
					
					#-- it is no necessary to send the manelist dictionary, just send the out_dir path
					img_save_dir, mod_save_dir = FILDA_IO.ini_output_dir(namelist, time)
					version = 'uiowa003'
					savename = img_save_dir + sat + 'IMG' + '.' + jdn + '.' +  overpass + '.' + version +  '.nc'
					# savename = './OUTPUT/' + sat + '.' + 'IMG' + '.' + jdn + '.' +  overpass + '.' + version +  '.nc'
					FILDA_IO.write_nc(fire_img, sat, savename)

					savename = mod_save_dir + sat + 'MOD' + '.' + jdn + '.' +  overpass + '.' + version +  '.nc'
					# savename = './OUTPUT/' + sat + '.' + 'MOD' + '.' + jdn + '.' +  overpass + '.' + version +  '.nc'
					FILDA_IO.write_nc(fire_mod, sat, savename)

else:
	print('The process could not start due to the required input files are not complete.')
