"""
Function to download the GEOS-FP from GMAO


import GEOS_FP
import os

base_dir = './'
date_work = "2022-09-01"

# retrieve geos url files for a given day
geos_urls = GEOS_FP.tavg1_url(date_work.replace('-',''), collections = ['tavg1_2d_flx_Nx'])

urls = []

# orbits to process 
overpass_to_process = ['0400']

# search on the urls file the one matching with the overpass
for file_url in geos_urls :
	geos_time = os.path.basename(file_url).split('_')[4][0:2]

	flag_time = [ True for orb_hr in overpass_to_process if orb_hr[0:2] == geos_time]
	if flag_time:
		urls.append(file_url)

# download the files matching with the overpass			
for url in urls:
	print(urls)
	bname_file = os.path.basename(url)
	command = 'wget ' + url + ' -O ' + base_dir +'/'+ bname_file
	os.system(command + "> /dev/null 2>&1")

"""

import datetime
import os
import subprocess

#------------------------------------------------------------------------------
#
def tavg1_url(yyyymmdd, collections=None):
    """ Generate a list that contains the URLs of GEOS-FP 
    one-day asm tavg1 data if collections is None.
    But users can specific variables now.
    (ywang, 03/17/2020)

    Parameters
    ----------
    yyyymmdd : str
        Date, '20180501' for example
    collections : list
        list of variables.  

    Returns
    -------
    urls : list
        All URLs
    """

    # URL of GEOS-FP data
    url = 'https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/'

    #
    if collections is None:
        collections = ['tavg1_2d_flx_Nx', 'tavg1_2d_lnd_Nx',
                       'tavg1_2d_rad_Nx', 'tavg1_2d_slv_Nx']

    # Generate URLs
    urls = []
    yyyy = yyyymmdd[0:4]
    mm   = yyyymmdd[4:6]
    dd   = yyyymmdd[6:8]
    url_curr_day = url + 'Y' + yyyy + '/M' + mm + '/D' + dd + '/'

    # tavg1
    Nh = 24
    for ih in range(Nh):

        ch = str(ih).zfill(2) + '30'

        for coln in collections:

            if 'tavg1' in coln:

                filename = url + 'Y' + yyyy + '/M' + mm + '/D' + dd + \
                        '/GEOS.fp.asm.' + coln  + '.' + yyyymmdd  + '_' + \
                        ch + '.V01.nc4'

                urls.append(filename)

    # inst3
    Nh = 8
    for ih in range(Nh):

        ch = str(ih*3).zfill(2) + '00'

        for coln in collections:

            if 'inst3' in coln:

                filename = url + 'Y' + yyyy + '/M' + mm + '/D' + dd + \
                        '/GEOS.fp.asm.' + coln  + '.' + yyyymmdd  + '_' + \
                        ch + '.V01.nc4'

                urls.append(filename)

    return urls

