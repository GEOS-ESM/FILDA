'''
This is the small PGE for generating the 90 NTL climatology


python FILDA_Gen_NTL_main.py 2019 08 h10v05

The three arguments are year, month, and tile 


'''

import glob
import numpy as np
import os
import sys
from datetime import datetime

import FILDA_IO			# module for reading and writing files
import FILDA_Time_Cord  # module for manipulating the time and coordinates
import FILDA_Resample 	# modele for physical resampling from DNB to M/I band
import FILDA_Cloud 		# modele for cloud masking
import FILDA_CLT	    # modele for generating the surface light climatology
import FILDA_NTL


year    = sys.argv[1]
month 	= sys.argv[2]
tile     = sys.argv[3]

# year  = '2019'
# month = '08' 
# tile  = 'h10v05'

data_Dir = '/Dedicated/jwang-data/shared_satData/VIIRS/yearly_data/VNP46A1/'
output_dir = '../Static/NTL/'

prefix = 'VNPNTL'
version = '001'

jdn = FILDA_Time_Cord.JulianDay(year + month+'01', outtype = 'nasa')

now = datetime.now()
current_time = now.strftime("%Y%m%d%H%M%S")

# According to Carol the name convention should be
# VNPNLT.AYYYYDDD.hxxvyy.vvv.YYYYMMDDHHMMSS.nc
savename = prefix + '.' + jdn + '.' + tile + '.' + version + '.' + current_time + '.nc'

# call the function to generate the climatology
FILDA_NTL.gen_NTL_climatology( year + month, tile, data_Dir, output_dir, savename)