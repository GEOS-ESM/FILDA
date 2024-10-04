'''

Python library of Multichannel Biphasic Fire Parameter (BPFPE)

space and time module, v1.0

This module provide the basic functions for coordinates and time 
manipulating.

'''


# function of date conversion
#-----------------------------------------------------------------------
def JulianDay(GregDay, **kwargs):
	'''
	Function to convert the Gregorian date to Julian Day Number(JDN)
	
	Reference: https://en.wikipedia.org/wiki/Julian_day

	Parameters
	----------
	GregDay	: string
			  format - yyyymmdd
			  
	Optional Parameters
	----------
	outtype : string, specify the output type, int, string, or nasa format
	type	: string
			  global - output is global Julian Day Number since November 23, −4713
			  local  - output is localized to that year
	Return
	----------    
	jdn		: int or str, defalt type int
			  Julian Day Number
	'''
	import numpy as np
	outtype = kwargs.get('outtype', 'int')
	type = kwargs.get('type', 'global')


	year = int(GregDay[0:4])
	month = int(GregDay[4:6])
	day = int(GregDay[6:8])
	Julian_a = (14-month)//12
	Julian_y = year + 4800 - Julian_a
	Julian_m = month + 12 * Julian_a - 3


	jdn = day + (153*Julian_m+2)//5 + 365*Julian_y + Julian_y//4 - \
	      Julian_y//100 + Julian_y//400 -32045


	if type == 'local':
		jdn_base = JulianDay( str(year) + '0101', type = 'global' )
		jdn = jdn - jdn_base + 1

	if outtype == 'str':
		jdn = np.str(jdn)

	if outtype == 'nasa':
		jdn_base = JulianDay( str(year) + '0101', type = 'global' )
		jdn = jdn - jdn_base + 1
		jdn = '00000' + np.str(jdn)
		jdn = 'A' + GregDay[0:4] + np.str(jdn)[-3:]
	
	
	return jdn

#-----------------------------------------------------------------------
def GregorianDay(jnd, outputformat = 'yyyy-mm-dd', **kwargs):
    '''
    Function to convert the Julian Day Number(JDN)the Gregorian date
    http://aa.usno.navy.mil/faq/docs/JD_Formula.php

	Parameters
	----------
	jdn				: julian day of year
	outputformat	: number of day of the year

	Return
	----------    
	GregDay			: Gregorian date
	
    '''
    
    noYear = kwargs.get('noYear', False)
    l = jnd + 68569
    n = 4*l//146097
	
    l= l-(146097*n+3)//4
    year= 4000*(l+1)//1461001
    l= l-1461*year//4+31
    month= 80*l//2447
    day = l-2447*month//80
    l= month//11
    month= month + 2 - 12*l
    year= 100*(n-49) + year + l

    if day < 10:
        day = '0' + str(day)
    if month <10:
        month = '0' + str(month)
    GregDay = outputformat
    GregDay = GregDay.replace('dd',str(day))
    GregDay = GregDay.replace('mm',str(month))
    if noYear:
    	GregDay = GregDay.replace('yyyy-','')
    else:
    	GregDay = GregDay.replace('yyyy',str(year))

    return GregDay

#-----------------------------------------------------------------------
def get_displatDate(date, format = 'mm-dd-yyyy', **kwargs):
	
	
	'''
	
	get_displatDate converts the julian day number string into normal date string
	
	Parameters:
	----------
	date			: julian day of year, A2020135
	format			: number of day of the year

	Return:
	----------    
	displatDate		: normal date string
	
	
	'''
	
	year = date[1:5]
	yearlyJDN = date[5:]	
	jdbBase = year + '0101'
	globalJDN = JulianDay(jdbBase,outtype = 'int') + int(yearlyJDN) - 1
	displatDate = GregorianDay(globalJDN, outputformat = format, **kwargs)
	
	return displatDate


#-----------------------------------------------------------------------
def get_date_series(gregDayBeg, gregDayEnd, **kwargs):

	import numpy as np
	outtype = kwargs.get('outtype', 'jdn')
	outputformat = kwargs.get('outputformat', 'yyyy-mm-dd')
	
	jdnBeg = JulianDay(gregDayBeg)
	jdnEnd = JulianDay(gregDayEnd) + 1
	jdns = np.arange(jdnBeg, jdnEnd)
	
	dateSeries = []
	if outtype == 'jdn':
		dateSeries = list( jdns )
		
	if outtype == 'greg':
		for jdn in jdns:
			GregDay = GregorianDay(jdn, outputformat = outputformat)
			dateSeries.append(GregDay)
			
	if outtype == 'nasa':
		for jdn in jdns:
			GregDay = GregorianDay(jdn, outputformat = 'yyyymmdd')
			dateSeries.append(JulianDay(GregDay, outtype = 'nasa'))		
	
	return dateSeries

#-----------------------------------------------------------------------
class init_time():
	'''
	init_time defines the time class for MCBEF 
	'''
	
	def __init__(self, jdn, overpass_beg, overpass_end):
		self.GD 	= get_displatDate(jdn)
		self.JDN	= jdn
		self.DOY	= jdn[5:]
		self.OP1	= overpass_beg
		self.OP2	= overpass_end
		self.Y		= self.GD[6:]
		self.M		= self.GD[0:2]
		self.D		= self.GD[3:5]
# 		self.h		= self.OP1[0:2]
# 		self.m		= self.OP1[2:]
		self.GEOS	= get_displatDate(jdn, format = 'yyyymmdd')
		
	def __str__(self):
		attrs = [f"{key}: {getattr(self, key)}" for key in 
				 sorted(self.__dict__.keys())]
		return "\n".join(attrs)

#-----------------------------------------------------------------------
def get_time_str(jdn, overpass):
	
	time_dict = {}
	display_time = get_displatDate(jdn)
	year     = display_time[6:]
	month    = display_time[0:2]
	day      = display_time[3:5]
	hour 	 = overpass[0:2]
	mins 	 = overpass[2:]
	
	# get the time string for GEOS-FP..
	gregDay = get_displatDate(jdn, format = 'yyyymmdd')
	geos_time = gregDay + '_' + hour
	
	time_dict['year']      = year
	time_dict['month']     = month
	time_dict['day'] 	   = day
	time_dict['geos_time'] = geos_time
	time_dict['greg_day']  = geos_time
	
	return time_dict

# function of coordinates manipulating
# for PlateCarree (Linear lat/lon grid)
#-----------------------------------------------------------------------
def get_point_tile_PlateCarree(cord, pos):
	'''
	
	Function to calculate the tile names for a given point coordinate

	'''
	import numpy as np
	lat = cord[0]
	lon = cord[1]

	vid = 8 - lat // 10
	hid = 18 + lon // 10
	
	res_lat = lat%10
	res_lon = lon%10
	
	if res_lat == 0:
		if pos == 'UpperLeft':
			vid = np.int(vid) + 1
		if pos == 'LowerRight':
			vid = int(vid)
		if pos == 'UpperRight':
			vid = int(vid) + 1
		if pos == 'LowerLeft':
			vid = int(vid)

	if res_lon == 0:
		if pos == 'UpperLeft':
			hid = int(hid)
		if pos == 'LowerRight':
			hid = int(hid) - 1
		if pos == 'UpperRight':
			hid = int(hid) - 1	
		if pos == 'LowerLeft':
			hid = int(hid)
	
	hid = int(hid)
	vid = int(vid)

	strhid = str(hid)
	
	while len(strhid) < 2:
		strhid = '0' + strhid

	strvid = str(vid)
	while len(strvid) < 2:
		strvid = '0' + strvid
		
	tile = 'h' + strhid + 'v' + strvid
	
	return hid, vid, tile

#-----------------------------------------------------------------------
def get_tiles(cord):
	
	'''
	
	Function to calculate the tile names for the given coordinates

	'''	
	
	import numpy as np
	
	hid_top, vid_top, _ = get_point_tile_PlateCarree((cord[0], cord[2]), pos = 'UpperLeft')
	hid_bot, vid_bot, _ = get_point_tile_PlateCarree((cord[1], cord[3]), pos = 'LowerRight')

	hids = np.arange(hid_top, hid_bot + 1, 1)
	vids = np.arange(vid_top, vid_bot + 1, 1)

	tiles = []
	for hid in hids:
		for vid in vids:
			strhid = str(hid)
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(vid)
			while len(strvid) < 2:
				strvid = '0' + strvid
		
			tile = 'h' + strhid + 'v' + strvid
			tiles.append(tile)
	return tiles

#-----------------------------------------------------------------------
def get_cord_PlateCarree(cord, numCeil = 2400, debug = 0):
	'''
	Function to get the geographical coordinates of the given region. 
	The coordinates are correspoding to Sinusodial grid.
	
	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates 
			   (Top Latitude, Bottom Latitude, Left Longitude, Right Longitude)
		
		numCeil : number of pixel in one tile, default value 1200
		
	Return
    ----------
		meshLat : array
				  longitude 
		meshLon : array
				  longitude 
		hidMin	: int
    		      vetical index of the point 
		hidMax	: int
    		  	 vetical index of the point 
		vidMin  : int
    		  	  vetical index of the point 
		vidMax  : int
    		  	  vetical index of the point
	'''
	import numpy as np

	UpperLeft = (cord[0], cord[2])
	UpperRight = (cord[0], cord[3])
	LowerLeft = (cord[1], cord[2])
	LowerRight = (cord[1], cord[3])

	tileInfor = []

		
	tileInfor.append(get_point_tile_PlateCarree(UpperLeft, 'UpperLeft'))
	tileInfor.append(get_point_tile_PlateCarree(UpperRight, 'UpperRight'))
	tileInfor.append(get_point_tile_PlateCarree(LowerRight, 'LowerRight'))
	tileInfor.append(get_point_tile_PlateCarree(LowerLeft, 'LowerLeft'))

	hid = []
	vid = []
	for item in tileInfor:
		hid.append(item[0])
		vid.append(item[1])
	
	hidMax = np.max(hid)
	hidMin = np.min(hid)

	vidMax = np.max(vid)
	vidMin = np.min(vid)
	
	
	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1

	# update the hid and vid
	hid = np.arange(hidMin, hidMax + 1, 1)
	vid = np.arange(vidMin, vidMax + 1, 1)
	
	GridDim = (num_v * numCeil, num_h * numCeil)
	meshLat = np.full(GridDim, np.nan)
	meshLon = np.full(GridDim, np.nan)
	
	if debug == 1:
		print('\n  - get_cord_VNP46A1 - hid: ', hid)
		print('\n  - get_cord_VNP46A1 - vid: ', vid)
	

	for hh in hid:
		for vv in vid:	
			strhid = str(hh)
			while len(strhid) < 2:
			
				strhid = '0' + strhid

			strvid = str(vv)
			while len(strvid) < 2:
				strvid = '0' + strvid
	
			tile = 'h' + strhid + 'v' + strvid
	
			if debug == 1:
				print('\n  - get_cord_VNP46A1 - tile: ',tile, hh, vv)
			
			latitude, longitude = cal_PlateCarree_grid(tile, numCeil)
			
			
			if debug == 1:
				print('\n  - get_cord_VNP46A1 - latitude: ',latitude[0,0], latitude[-1,0])
			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			
			if debug == 1:
				print(hIdx, vIdx)
				print(vIdx * numCeil, (vIdx + 1) * numCeil, hIdx * numCeil, (hIdx + 1) * numCeil)
		
			meshLat[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = latitude
				
			meshLon[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = longitude	
					
	return meshLat, meshLon, hidMin, hidMax, vidMin, vidMax

#-----------------------------------------------------------------------
def cal_PlateCarree_grid(tile, numCeil):
	
	'''
	Function of getting the coordinates of the Plate-Carree grid
	
	
	'''
	
	import numpy as np
	
	vid = np.int(np.float(tile[4:6]))
	hid = np.int(np.float(tile[1:3]))

	latBoundary = [(8 - vid) * 10, (9 - vid) * 10]
	lonBoundary = [(hid - 18) * 10, (hid - 17) * 10]
	
	latitude  = np.linspace(latBoundary[1], latBoundary[0], numCeil)

	longitude  = np.linspace(lonBoundary[0], lonBoundary[1], numCeil)

	latitude = (latitude * np.ones((numCeil,1),np.float32)).T
	
	longitude = np.ones((numCeil,1),np.float32) * longitude
    
	return latitude, longitude


# for Sinusoidal (Linear lat/lon grid)
#-----------------------------------------------------------------------
def cal_sinu_xy(tile, numCeil):
    '''

    Function to calculate the geographical coordinates of the sinusoidal grid

    Parameters
    ----------
        tile - str format, example: 'h07v05'
        numCeil - number of the ceils in one tile

    Return
    ---------- 
        geographical coordinates of the tile

    Reference: 1. https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
               2. https://onlinelibrary.wiley.com/doi/pdf/10.1111/0033-0124.00327
               3. https://modis-land.gsfc.nasa.gov/GCTP.html
    
    MODIS use 6371007.181 as the radius of the Earth...
    '''
    import numpy as np
    
    numHoriTail = 37
    numVertTail = 19
    
    halfHoriLenght = geog_to_sinu([0, 180])[0]
    halfVertLenght = geog_to_sinu([90, 0])[1]
    resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
    halfCeilLen = resol_ceil/2.0



    xx = np.linspace(-halfHoriLenght, halfHoriLenght, numHoriTail)
    yy = np.linspace(halfVertLenght, -halfVertLenght, numVertTail) 

    vid = int(float(tile[4:6]))
    hid = int(float(tile[1:3]))
    #print('  - Calculating sinusoidal coordinates of ', tile)
    #print('    Vertical Tile:', vid, 'Horizontal Tile:', hid)

    x = np.linspace(xx[hid], xx[hid+1], numCeil + 1)
    y = np.linspace(yy[vid], yy[vid+1], numCeil + 1)
    
    x = (x[0:-1] + x[1:])/2.
    y = (y[0:-1] + y[1:])/2.
    
    xv, yv = np.meshgrid(x, y)
    
    return xv, yv, resol_ceil

#-----------------------------------------------------------------------
def cal_sinu_grid(tile, numCeil):
    '''

    Function to calculate the geographical coordinates of the sinusoidal grid

    Parameters
    ----------
        tile - str format, example: 'h07v05'
        numCeil - number of the ceils in one tile

    Return
    ---------- 
        geographical coordinates of the tile

    Reference: 1. https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
               2. https://onlinelibrary.wiley.com/doi/pdf/10.1111/0033-0124.00327
               3. https://modis-land.gsfc.nasa.gov/GCTP.html
    
    MODIS use 6371007.181 as the radius of the Earth...
    '''
    
    import numpy as np
    
    numHoriTail = 37
    numVertTail = 19
    
    halfHoriLenght = geog_to_sinu([0, 180])[0]
    halfVertLenght = geog_to_sinu([90, 0])[1]
    resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
    halfCeilLen = resol_ceil/2.0



    xx = np.linspace(-halfHoriLenght, halfHoriLenght, numHoriTail)
    yy = np.linspace(halfVertLenght, -halfVertLenght, numVertTail) 

    vid = int(float(tile[4:6]))
    hid = int(float(tile[1:3]))
    #print('  - Calculating sinusoidal coordinates of ', tile)
    #print('    Vertical Tile:', vid, 'Horizontal Tile:', hid)

    x = np.linspace(xx[hid], xx[hid+1], numCeil + 1)
    y = np.linspace(yy[vid], yy[vid+1], numCeil + 1)
    
    x = (x[0:-1] + x[1:])/2.
    y = (y[0:-1] + y[1:])/2.
    
    xv, yv = np.meshgrid(x, y)
    
    latitude, longitude = sinu_to_geog((xv, yv))

    return latitude, longitude, resol_ceil


#-----------------------------------------------------------------------
def sinu_to_geog(cord):
	'''
	Function of converting the sinusoidal projection to platecree projection

	Parameters
	----------
		sinusoidal point - list or tuple liked, (x, y)

	Return
	----------
		geographical coordinates - latitude, longituede
	'''
	import numpy as np
	x = cord[0]
	y = cord[1]
	pi = 180.0 / np.pi
	R = 6371007.181000

	phi = y/R
	lamda = x / np.cos(phi) / R

	latitude = phi * pi
	longituede = lamda * pi

	return latitude, longituede

#-----------------------------------------------------------------------
def geog_to_sinu(cord):
    '''
    Function of converting the geographical projection to sinusoidal projection

    Parameters
    ----------
        geographical coordinates - list or tuple liked, (latitude, longituede)

    Return
    ----------
        sinusoidal point - (x, y)


    '''

    import numpy as np

    lat = cord[0]
    lon = cord[1]

    pi = 180.0 / np.pi
    R = 6371007.181000

    phi = lat / pi
    lamda = lon / pi
    y = phi * R
    x = np.cos(phi) * lamda * R

    return x, y

#-----------------------------------------------------------------------
def get_tile_sinusoidal(cord, namelist):


	import numpy as np
	
	#select sinusoidal tile according the txt
	sinu_data=np.loadtxt(namelist['LUT_DIR'] + 'sinusoidal_tile.txt')
	lon_min=sinu_data[:,2]
	lon_max=sinu_data[:,3]
	lat_min=sinu_data[:,4]
	lat_max=sinu_data[:,5]    

	UpperLeft = (cord[0], cord[2])
	UpperRight = (cord[0], cord[3])
	LowerLeft = (cord[1], cord[2])
	LowerRight = (cord[1], cord[3])

	corners = [UpperLeft, UpperRight, LowerLeft, LowerRight]

	hid = []
	vid = []
	for C in corners:
		i = 0
		in_tile = False
		while(not in_tile):
			in_tile = C[0] >= sinu_data[i, 4] and C[0] <= sinu_data[i, 5] and C[1] >= sinu_data[i, 2] and C[1] <= sinu_data[i, 3]
			i += 1
		vert = sinu_data[i-1, 0]
		horiz = sinu_data[i-1, 1]
	
		hid.append(horiz)
		vid.append(vert)
	
		print( ' - ', C[0], C[1],' Vertical Tile:', vert, 'Horizontal Tile:', horiz)
	hidMax = int(np.max(hid))
	hidMin = int(np.min(hid))

	vidMax = int(np.max(vid))
	vidMin = int(np.min(vid))
	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1

	# update the hid and vid
	hids = np.arange(hidMin, hidMax + 1, 1)
	vids = np.arange(vidMin, vidMax + 1, 1)    

	tiles = []
	for hid in hids:
		for vid in vids:
			strhid = str(int(hid))
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(int(vid))
			while len(strvid) < 2:
				strvid = '0' + strvid

			tile = 'h' + strhid + 'v' + strvid
			tiles.append(tile)

	return tiles, hidMax, hidMin, vidMax, vidMin

#-----------------------------------------------------------------------
def get_point_in_tile_Sinusoidal(cord, numCeil = 2400):
	'''
	Function of calculating the tile of a specific point

	Parameters
	----------
		cord : list or tuple
			   geographical coordinates (latitude, longituede)

		pos : position of the point

	Return
	----------
		hid : int
			  horizental index of the point
		vid : int
			  vetical index of the point
		tile: str
			  tile name
	'''
	import numpy as np


	numHoriTail = 37
	numVertTail = 19

	halfHoriLenght = geog_to_sinu([0, 180])[0]
	halfVertLenght = geog_to_sinu([90, 0])[1]

	tileHoriLenght = halfHoriLenght / (numHoriTail-1) * 2
	tileVertLenght = halfVertLenght / (numVertTail-1) * 2

	resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
	halfCeilLen = resol_ceil/2.0    

	# process the four corner

	x, y = geog_to_sinu( cord )


	x_res = abs(np.round(x / tileHoriLenght) * tileHoriLenght - x)
	y_res = abs(np.round(y / tileVertLenght) * tileVertLenght - y)    
	hid =  np.round(x // tileHoriLenght) + 18
	vid = 8 -  np.round(y // tileVertLenght)


	return int(hid), int(vid)
    
#-----------------------------------------------------------------------
def get_tile_sinusoidal_2(cord, namelist, numCeil = 2400):
    
	import numpy as np
	#select sinusoidal tile according the txt

	numHoriTail = 37
	numVertTail = 19

	halfHoriLenght = geog_to_sinu([0, 180])[0]
	halfVertLenght = geog_to_sinu([90, 0])[1]

	tileHoriLenght = halfHoriLenght / (numHoriTail-1) * 2
	tileVertLenght = halfVertLenght / (numVertTail-1) * 2

	resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
	halfCeilLen = resol_ceil/2.0    

	tileInfor = {}
	tileInfor['UpperLeft'] = (cord[0], cord[2])
	tileInfor['UpperRight'] = (cord[0], cord[3])
	tileInfor['LowerLeft'] = (cord[1], cord[2])
	tileInfor['LowerRight'] = (cord[1], cord[3])

	# process the four corner  
	hids = []
	vids = []

	#----
	hid, vid = get_point_in_tile_Sinusoidal(tileInfor['UpperLeft'], numCeil = numCeil)
	hid = hid - 1
	hids.append(hid)  
	vid = vid - 1
	vids.append(vid)

	#----
	hid, vid = get_point_in_tile_Sinusoidal(tileInfor['UpperRight'], numCeil = numCeil)
	hid = hid + 1
	hids.append(hid)   
	vid = vid - 1
	vids.append(vid)   

	#----
	hid, vid = get_point_in_tile_Sinusoidal(tileInfor['LowerLeft'], numCeil = numCeil)   
	hid = hid - 1
	hids.append(hid) 
	vid = vid + 1
	vids.append(vid)  

	#----
	hid, vid = get_point_in_tile_Sinusoidal(tileInfor['LowerRight'], numCeil = numCeil)
	hid = hid + 1
	hids.append(hid)  
	vid = vid + 1
	vids.append(vid)      


	hidMax = int(np.max(hids))
	hidMin = int(np.min(hids))

	vidMax = int(np.max(vids))
	vidMin = int(np.min(vids))


	if hidMin <=0:
		hidMin = 0
	if hidMax >=35:
		hidMax = 35
	

	if vidMin <=0:
		vidMin = 0
	if vidMax >=17:
		vidMax = 17

	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1

	# update the hid and vid
	hids = np.arange(hidMin, hidMax + 1, 1)
	vids = np.arange(vidMin, vidMax + 1, 1)

	tiles = []
	for hid in hids:
		for vid in vids:
			strhid = str(int(hid))
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(int(vid))
			while len(strvid) < 2:
				strvid = '0' + strvid

			tile = 'h' + strhid + 'v' + strvid
			tiles.append(tile)

	return tiles, hidMax, hidMin, vidMax, vidMin

#-----------------------------------------------------------------------
def get_point_in_tile_Sinusoidal_3(cord, numCeil = 2400):
	'''
	Function of calculating the tile of a specific point

	Parameters
	----------
		cord : list or tuple
			   geographical coordinates (latitude, longituede)

		pos : position of the point

	Return
	----------
		hid : int
			  horizental index of the point
		vid : int
			  vetical index of the point
		tile: str
			  tile name
	'''
	import numpy as np


	numHoriTail = 37
	numVertTail = 19

	halfHoriLenght = geog_to_sinu([0, 180])[0]
	halfVertLenght = geog_to_sinu([90, 0])[1]

	tileHoriLenght = halfHoriLenght / (numHoriTail-1) * 2
	tileVertLenght = halfVertLenght / (numVertTail-1) * 2

	resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
	halfCeilLen = resol_ceil/2.0    

	# process the four corner

	x, y = geog_to_sinu( cord )


	x_res = abs(np.round(x / tileHoriLenght) * tileHoriLenght - x)
	y_res = abs(np.round(y / tileVertLenght) * tileVertLenght - y)    
	hid =  np.round(x // tileHoriLenght) + 18
	vid = 8 -  np.round(y // tileVertLenght)


	return hid.astype(int), vid.astype(int)


#-----------------------------------------------------------------------
def get_tile_sinusoidal_3(cord, numCeil = 2400):

	import numpy as np
	#select sinusoidal tile according the txt

	numHoriTail = 37
	numVertTail = 19

	halfHoriLenght = geog_to_sinu([0, 180])[0]
	halfVertLenght = geog_to_sinu([90, 0])[1]

	tileHoriLenght = halfHoriLenght / (numHoriTail-1) * 2
	tileVertLenght = halfVertLenght / (numVertTail-1) * 2

	resol_ceil = halfHoriLenght/((numHoriTail-1)/2.)/numCeil
	halfCeilLen = resol_ceil/2.0    


	hids, vids = get_point_in_tile_Sinusoidal_3(cord, numCeil = numCeil)

	
	hids = list(set(hids.tolist()))
	vids = list(set(vids.tolist()))
	
	hidMax = np.max(hids) + 1
	hidMin = np.min(hids) - 1

	vidMax = np.max(vids) + 1
	vidMin = np.min(vids) - 1

	if hidMin <=0:
		hidMin = 0
	if hidMax >=35:
		hidMax = 35
	

	if vidMin <=0:
		vidMin = 0
	if vidMax >=17:
		vidMax = 17

	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1
	

	# update the hid and vid
	hids = np.arange(hidMin, hidMax + 1, 1)
	vids = np.arange(vidMin, vidMax + 1, 1)
	
	tiles = []
	for hid in hids:
		for vid in vids:
			strhid = str(int(hid))
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(int(vid))
			while len(strvid) < 2:
				strvid = '0' + strvid

			tile = 'h' + strhid + 'v' + strvid
			tiles.append(tile)

	
	return tiles, hidMax, hidMin, vidMax, vidMin

#-----------------------------------------------------------------------
# Function to convert HHMM string to the i-th 6-minute interval index
def convert_to_interval_index(times):
    # Extract hours and minutes from the HHMM format
    hours = times.astype(int) // 100   # Integer division by 100 to get hours
    minutes = times.astype(int) % 100  # Modulo 100 to get minutes

    # Calculate total minutes from the start of the day
    total_minutes = hours * 60 + minutes

    # Calculate the index for 6-minute intervals, starting count from 1
    interval_index = (total_minutes // 6)

    return interval_index

