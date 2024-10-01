# FILDA2
The Fire Light Detection Algorithm Version 2 (FILDA-2)

VNP47/VJ147 is a new experimental Level-2 (L2) VIIRS nighttime active fire combustion efficiency products at both M-band (750m) and I-band (375m) resolution. The product’s algorithm is based on the most recent two papers on nighttime fire studies (Polivka et al., 2016, Wang et al., 2020). In essence, it has its own cloud screening, water body removal process, and an approach to integrate the I-band, M-band and DNB for fire characterization. The products provide the fire location (latitude/longitude), fire radiative power (FRP), visible energy fraction (VEF), and modify combustion efficiency (MCE) of  fire at the pixel level globally, every day.​

##  - **Recent update** - 

1. Fix rounding problem when calculating the pixel index of the detected fire

## Table of content

- [Package Dependence](#package-dependence)

- [How to use FILDA-2](#how-to-use-filda-2)

- [Namelist](#namelist)

- [Input data](#input-data)
  - [Daily updated data](#daily-updated-data)
  - [Three month updated data](#three-month-updated-data)
  - [Yearly updated data](#yearly-updated-data)
  - [Static data](#static-data)

- [Layout of the input dataset](#layout-of-the-input-dataset)

- [Output](#output)

- [Layout of the output dataset](#layout-of-the-output-dataset)

- [Package workflow](#package-workflow)

- [Create the nighttime climatology](#create-the-nighttime-climatology)

- [Download GEOS-FP](#download-geos-fp)

- [Runtime](#runtime)

- [Institution](#institution)

- [Contact](#contact)

- [Website](#website)



## Package Dependence
 The FILDA-2 python package is developed on python (3.7.9)
 FILDA-2 dependents on the following python libraries:

 
| #  |    libraries     | version   | Description| 
| ---|       :---:      |   :---:   | ------------- |
| 1  | os               | Default   | Miscellaneous operating system interfaces |
| 2  | sys              | Default   | System-specific parameters and functions  |
| 3  | copy             | Default   | Shallow and deep copy operations  |
| 4  | datetime         | Default   | Basic date and time types  |
| 5  | multiprocessing  | Default   | Process-based parallelism  |
| 6  | numpy            | 1.19.1    | https://numpy.org/  |
| 7  | scipy            | 1.4.1     | https://scipy.org/  |
| 8  | pandas           | 1.1.3     | https://pandas.pydata.org/	  |
| 9  | netCDF4          | 1.5.1.2   | https://unidata.github.io/netcdf4-python/  |
| 10 | h5py             | 3.7.0     | https://docs.h5py.org/en/stable/build.html#source-installation |
| 11 | pyhdf            | 0.10.5    | https://hdfeos.org/software/pyhdf.php |

```FILDA.yml``` contains the list of the necessary libraries to run FILDA-2. We recommend to create an environment and
install all of these libraries at the time. The above can be performed by using ```Anaconda``` or ```Pip``` as follows.

### Anaconda
If you work under Anaconda (preferable, a library manage tool that help to check the interdenpendency of libraries), you can create an anaconda 
environment to install all the libraries by using the following anaconda command:
```
conda env create -f FILDA.yml
```

To run FILDA-2 in the FILDA python environment, you need to activate the environment using the following command:

```
source activate FILDA
```
### Pip
Alternatively, the list of FILDA-2 requirements can be installed into a standard python virtual environment and then 
installing the libraries by using the ```pip``` packages insalled. The general steps for this approach are:
1. Create a new environment
```
    python -m venv env
```
2. Activate the virtual environment
```
    source env/bin/activate
```
3. Install the libraries specified in the FILDA.yml file
```
    pip install -r FILDA.yml
```

In both cases, once the environment is created and the requirement are installed, you should be able to implement 
the FILDA-2 library.


## How to use FILDA-2:
FILDA-2 operates at overpass level, as input parameters they are required:
```
- jdn, julian day in the format "AYYYYDDD", e.g. "A2022006"
- overpass, time in the format "HHMM", e.g. "0018"
- sat, satellite, options: VNP | VJ1
```
Note: the parameters comes from the input file names, as an example see the below figure:
![alt text](https://github.com/menkbard/FILDA2/blob/master/docs/parameters.png)


For case of using a *bash script* to control input parameters for running the package, use script in the following manners:
```
FILDA_main.py A2022006 0018 VNP
```

For case of using a *python script* to control input parameters for running the package, use the script:
```  
FILDA_main_op.py
```

In the file ```test_run.py ``` file it is provided an example of FILDA-2 implementation

## Namelist
The namelist is the configuration file to control the FILDA2 package. There are three blocks in the current namelist, **Static dataset**, **Run & Output**, **Thresholds**. A operational user will only need to mofidy the first and second section. 

```
# FILDA2 v1.0 namelist file Oct. 10 2022, Meng Zhou
#-----------------------------------------------------------------------
# Static dataset
NTL_DIR             : ...
...
#-----------------------------------------------------------------------
# Run & Output
# User need to put the VIIRS Level1b and GEOS-FP data under the RUN_DIR
RUN_DIR             : ...
...
#-----------------------------------------------------------------------
# Adjustable thresholds
# Users do not change for operational purpose
thres_BTI04         : ...
```

### Static dataset
```
# FILDA2 v1.0 namelist file Oct. 10 2022, Meng Zhou
#-----------------------------------------------------------------------
# Static dataset
NTL_DIR             : Path of the VIIRS Nighttime Light (NTL) Climatology
SURF_DIR            : Path of the Land surface type dataset
GASFLARE_DIR        : Path of the annual gas flaring dataset
PEATLAND_DIR        : Path of the peatland dataset
LUT_DIR             : Path of the static look-up table (LUT) for FILDA2
# The default year for the land surface type and the nighttime climatology
# if FILDA2 cannot find the the correct static data provide, it will try to 
# use the following years.
Default_year_SURF   : The latest year of the Land surface type dataset available
Default_year_NTL    : The latest year of the NTL dataset available
# prefix for static data
NTL_prefix          : Prefix for the NTL dataset
PEAT_prefix         : Prefix for the Peatland dataset
```

In preparing the static data, some specific layout of each static dataset is required, please refer to [Layout of the output dataset](#layout-of-the-output-dataset) for the preparation of the static data. 

**Default_year_SURF** is the years that the latest land surface type dataset available. This is a compromise that the MCD12Q1 has a relative long time latency. In the workflow, FILDA2 will try to find the correct year for the land surface type, e.g 2017 land surface type for 2017 fire detection. However, for near real time processing, it is highly possible that when conducting fire detection for 2022, the latest land surface type dataset available will only be 2020. In this senario, FILDA2 will use the 2020 dataset as a proxy. 

**Default_year_NTL** same reason as the **Default_year_SURF** but for the nighttime light climatology.

In the **FILDA2_test_demo** provided by the AER Lab FTP, the static data are provided in minimal just for processing the level-1b gradule provided, the setting of the **Static dataset** in the namelist is
```
# FILDA2 v1.0 namelist file Oct. 10 2022, Meng Zhou
#-----------------------------------------------------------------------
# Static dataset
NTL_DIR             : ../Static/NTL/
SURF_DIR            : ../Static/MCD12Q1/
GASFLARE_DIR        : ../Static/GASFLARE/
PEATLAND_DIR        : ../Static/PEATLAND/
LUT_DIR             : ../Static/LUT/
# The default year for the land surface type and the nighttime climatology
# if FILDA2 cannot find the the correct static data provide, it will try to 
# use the following years.
Default_year_SURF   : 2020
Default_year_NTL    : 2019
# prefix for static data
NTL_prefix          : VNPNTL
PEAT_prefix         : N/A
```

### Run & Output directories
```
#-----------------------------------------------------------------------
# Run & Output
# User need to put the VIIRS Level1b and GEOS-FP data under the RUN_DIR
RUN_DIR             : Path of the VIIRS level1B and the GEOS-FP data
OUT_DIR             : Path of the FILDA2 output
Layout_MODE         : 1 or 3
COPY_STATIC         : 0 or 1
```

The **RUN_DIR** is the directory that holds the VIIRS level-1B and GEOS-FP data. 

The **OUT_DIR** is the root directory for the FILDA2 output. By default, FILDA2 will create an **IMG** and an **MOD** subdirectory to store the I-band and M-band output seperately. Please refer to [Layout of the output dataset](#layout-of-the-output-dataset) for detailed information.

The **Layout_MODE** is a switch determining the way of FILDA2 to find the correct VIIRS level1-b and GEOS-FP data

 - 1 : The layout of the VIIRS level-1b follows the LAADS, GEOS-FP data are stored in the **GEOS-FP** subdirectory under the **RUN** subdirectory.
 - 3 : The level-1b and GEOS-FP data (7 in total) are directory stored in the the **RUN** subdirectory

For more layout information, please refer to [Layout of the input dataset](#layout-of-the-input-dataset). 

The **COPY_STATIC** is a switch determing whether copy the necessary static data into the RUN directory.

 - 0 : Do not copy. FILDA2 will read static data from the directorie provided in the **Static dataset**
 - 3 : Copy. FILDA2 will copy a subset of thoese static data into the **RUN** directory and read those data from the **RUN** directory.


The following setting applies to the senario that the user prepares the VIIRS level-1b and GEOS-FP data (7 in total) in a directory called **RUN**. When excuting FILDA2, those necessary static data will also be copied from their directories provided in the **Static dataset** in to the **RUN** directory.

```
#-----------------------------------------------------------------------
# Run & Output
# User need to put the VIIRS Level1b and GEOS-FP data under the RUN_DIR
RUN_DIR             : ./RUN/
OUT_DIR             : ./RUN/
Layout_MODE         : 3
COPY_STATIC         : 1
```



## Input data

FILDA relays on several inputs to detect fires, they are different in the time latency.

### Daily updated data
**1. VIIRS real-time level-1b data (For both radiance and geolocation).**

*Radiance data*
   
- VIIRS Level-1B calibrated DNB radiance data (VNP02DNB for Suomi-NPP VIIRS, VJ102DNB for NOAA20 VIIRS)
- VIIRS Level-1B calibrated MOD radiance data (VNP02MOD, VJ102MOD)
- VIIRS Level-1B calibrated IMG radiance data (VNP02IMG, VJ102IMG)

*Geolocation*
   
- VIIRS the DNB geolocation product (VNP03DNB, VJ103DNB)
- VIIRS the M-band geolocation product (VNP02MOD, VJ102MOD)
- VIIRS the I-band geolocation product (VNP02IMG, VJ102IMG)


The VIIRS level-1b data can be retrieved from LAADS data archive https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/. 5200 is the current collection used for the FILDA2.


**2. The surface flux diagnostics (tavg1_2d_flx_Nx)**

GEOS model data:
Goddard Earth Observing System (GEOS) forward processing (FP) is used to provide the land surface temperature (tavg1_2d_flx_Nx) for the cloud mask in the FILDA-2 M-band cloud masking. It undergoes regular updates to incoperate the latest development for the GOES system. Global Modeling and Assimilation Office (GMAO) will continue to produce the GOES-FP product. The resolution of the GEOS-IT is 0.25 X 0.3125.

The Goddard Earth Observing System (GEOS) for Instrument Team (IT), former GEOS Forward Processing for Instrument Teams (FP-IT), is recommended for the reprocessing of the FILDA-2 data. GEOS-IT uses a "semi-frozen" GEOS system to ensure long-term continuity and reproducibility. The resolution of the GEOS-IT is 0.5 X 0.625 degrees.

The FILDA-2 package read_Geos_FP(filename, GeosVar) is capable of reading the spatial resolution of the GEOS file provided in the **RUN** directory. 

The GEOS-FP tavg1_2d_flx_Nx can be obtained from https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/

GMAO has not released the GEOS-IT yet and GEOS-FP-IT could be used as a substitute for the reprocessing.


### Three-month updated data
**3. Nighttime light climatology**

The nigttime light climatology is generated upon the standard S-NPP VIIRS Black Marble Nighttime Light product (VNP46A1). It only needs to be updated at the begining of each month and will be valid for the entire month. For example, at July 1st, 2023, a global nighttime light climatologies will be generated for the July 2023. The data used for this climatology will be from March 23rd, 2023 to May 23rd, 2023. Here timeframe of a week is applied to compensate the possible latency of the VNP46A1 data. The format of the dataset is netCDF4.

### Yearly updated data
**4. Land cover database**

The land surface type dataset used in FILDA2 comes from MODIS Terra and Aqua combined land cover product (MCD12Q1). The format of the dataset is netCDF4.

**5. Gas flaring database**

The Gas flaring dataset is a rasterized dataset. The format of the dataset is netCDF4 (~2G).

### Static data
**6. Peat Land database**

The Peat land dataset is a rasterized dataset. The format of the dataset is netCDF4 (~2G).

**7. Look-up Tables**
The following LUTs used in the FILDA2.

- Infrared.csv : Temperature and view geometries dependent threshold for cloud mask
- VNP_DNB2IMG_Resampling_Lookup_Table.nc: DNB to IMG resampling LUT for S-NPP VIIRS
- VNP_DNB2MOD_Resampling_Lookup_Table.nc: DNB to MOD resampling LUT for S-NPP VIIRS
- VNP_MOD2DNB_Resampling_Lookup_Table.nc: MOD to DNB resampling LUT for S-NPP VIIRS
- VNP_MOD_Pixelareas_Lookup_Table.nc: MOD pixel area LUT for S-NPP VIIRS

- VJ1_DNB2IMG_Resampling_Lookup_Table.nc: DNB to IMG resampling LUT for NOAA20 VIIRS
- VJ1_DNB2MOD_Resampling_Lookup_Table.nc: DNB to MOD resampling LUT for NOAA20 VIIRS
- VJ1_MOD2DNB_Resampling_Lookup_Table.nc: MOD to DNB resampling LUT for NOAA20 VIIRS
- VJ1_MOD_Pixelareas_Lookup_Table.nc: MOD pixel area LUT for NOAA20 VIIRS
- 
**=====================================**

### IMPORTANT NOTE for dataset

**=====================================**

The FILDA 2 package adopts a flexible approach in selecting the required datasets, such as Nighttime Lights (VNPNTL), Land Surface (MCD12Q1), Peatland, and Gasflaring. These datasets are chosen based on availability but are not mandatory for the process. In essence, FILDA2 evaluates the necessary files for processing a VIIRS granule and attempts to locate them within the provided directory specified in the namelist. If FILDA2 successfully finds a file, it will read and utilize it accordingly. However, if a file is not found, it will raise a warning instead of breaking the code or triggering an error. This is particularly relevant for ocean surface data, as there may be cases where no dataset is available from the original provider.


## Layout of the input dataset

**1. Level-1b layout**

FILDA2 is developed using VIIRS data retrieved from LAADS data archive. For the layout of the level-1b input data. It follows the layout of the LAADS as follow: 

 - Root directory
    - Product
      - Year
        - day of year
          - *.nc file

For example, in LAADS, the collection 5200 VNP02DNB dataset has the following structure:

https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP02DNB/2022/006/

  - VNP02DNB
    - 2021
    - 2022
      - 005
      - 006
        - ...
        - VNP02DNB.A2022006.0018.002.2022269020936.nc
        - ...
      - ...
    - ...
  - VNP03DNB


**2. GEOS-FP layout**

Their is no specific layout requirement for the GOES-FP data, FILDA2 package will read GEOS-FP data from the _**GEOS_DIR**_ provided in the namelist.

**3. Nighttime light layout**

Their is no specific layout requirement for the Nighttime light data, FILDA2 package will read GEOS-FP data from the _**NTL_DIR**_ provided in the namelist.

**4. Land surface type data layout**

The Land surface type data used in the current FILDA2 package is MODIS MCD12Q1 which can be retrieved from the LAADS. It follows the LAADS layout.

**5. Gasflaring layout**

Their is no specific layout requirement for the Nighttime light data, FILDA2 package will read GEOS-FP data from the _**GASFLARE_DIR**_ provided in the namelist.


**6. Peatland layout**

Their is no specific layout requirement for the Nighttime light data, FILDA2 package will read GEOS-FP data from the _**PEATLAND_DIR**_ provided in the namelist.

**7. Lookup table layout**

Their is no specific layout requirement for the Nighttime light data, FILDA2 package will read GEOS-FP data from the _**LUT_DIR**_ provided in the namelist.

**======================**

**=== IMPORTANT NOTE ===**

**======================**

At the begin of the deteion, the following code will automatically match the level-1b and GEOS-FP input  for the FILDA2 with 

```
file_dict = FILDA_IO.get_files(sat, namelist, time) 
```

If a user would like to copy all the necessay level-1b and GEOS-FP files in to a **RUN** directory, then execute the FILDA2 package, the following mofidication can force the code to match the level-1b and GEOS-FP files in the same directory

```
Layout_MODE           : 3
```

## Output
For each overpass, FILDA2 will generate a I-band (375m) and M-band (750m) output. The current name convention at UIowa is 

**Sensor.Resolution.JulianDay.Overpass.Version.nc**

Below is an exampe of the outputs of August 1st, 2019, UTC0200

VNP.IMG.A2019213.0200.001.nc

VNP.MOD.A2019213.0200.001.nc

The user can modify the name convention at the end of FILDA_main.py,
```

version = '001'
savename = img_save_dir + sat + '.' + 'IMG' + '.' + jdn + '.' +  overpass + '.' + version +  '.nc'

```

The possible name is:

VNP47.A2019213.0200.001.time_of_process.nc [**M-band resolution product**]

VNP47IMG.A2019213.0200.001.time_of_process.nc [**I-band resolution product**]

## Layout of the output dataset

The package will generate a OUTPUT directory under **OUT_DIR** following the LADDAS data structure:

  - OUT_DIR
    - MOD
      - Year
        - Day of the year
          - *.nc
    - IMG
      - Year
        - Day of the year
          - *.nc

To modify the output structure, users need to modify the _**ini_output_dir**_ in the FILDA_IO.py module

```
#-----------------------------------------------------------------------
def ini_output_dir(namelist, time):

	'''
	ini_output_dir creates the out put path for a given time

	'''
```

## Package workflow
The figure shows blows is the work flow for Layout_MODE = 3, COPY_STATIC = 0. This is the current setting at UIowa, FILDA2 will read static data in the local directory instead of copying them to the **RUN** directory.
![alt text](https://github.com/menkbard/FILDA2/blob/master/docs/UIowa_work_flow.png)

The figure shows blows is the work flow for Layout_MODE = 3, COPY_STATIC = 1. In this setting, FILDA2 will first copy the static data to the **RUN** directory and read it from there.
![alt text](https://github.com/menkbard/FILDA2/blob/master/docs/work_flow_option_2.png)


## Create the nighttime climatology

The nighttime climatology is generated using the NASA Black Marble VNP46A1 product, available through the LAADS system at https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A1/. The stand-alone module **FILDA_NTL.py** is utilized for generating this climatology.

The script produces a 3-month (90-day) climatology for August 2019, specifically targeting tile h10v05. A data time lag of 7 days is applied, meaning the data utilized for generating the August 2019 climatology spans from April 26, 2019 (A2019116) to July 25, 2019 (A2019206).

For users unfamiliar with the linear latitude/longitude grid, reference is recommended to Figure 2 in the [tile_info](https://viirsland.gsfc.nasa.gov/PDF/BlackMarbleUserGuide_v1.2_20210421.pdf). It is crucial to prepare and store global nighttime climatologies at the start of each month for the corresponding month. The 'h' value ranges from 00 to 35, and the 'v' value from 00 to 17.

```
python FILDA_Gen_NTL_main.py 2019 08 h10v05
```

Note that the layout of the VNP46A1 follows the NASA LAADS layout as referenced above.
Use needs to goin to the FILDA_Gen_NTL_main.py to specificy the output directory.

```
output_dir = '../Static/NTL/'
```

### Rules
1. The NTL climatology is generated based on the NASA Black Marble VNP46A1 product.
2. The NTL climatology should be generated at the beginning of each month and remain valid for the entire month.
3. Users need to prepare the global nighttime climatologies at the start of each month and store them for that month. (h ranges from 00 to 35, v ranges from 00 to 17)
4. The climatology generated from VNP data can be applied to VJ1/VJ2 fire detection.
5. For January to April 2012, as there is no data preceding these months, the climatology for these months can be generated based on data from the months themselves (e.g., January 2012 to April 2012).
6. The gen_NTL_climatology function will handle the cross-year issue.


Function, **update_params_for_collection(filename, params)**, has been integrated into FILDA_NTL. This function scrutinizes the provided filename, ascertains the collection of the input VNP46A1, and modifies the **datafiled** and variable names for subsequent processes.

Function, **special_handling(year_month, tile, sat)**, has been introduced to FILDA_NLT. This function offers hardcoded solutions for generating climatology during the initial months of VNP and VJ1 and also addresses challenges in high latitude regions:
  
  - For VNP: The climatology from May 2012 is utilized to represent the climatology from January to April 2012.
  - For VJ1: The climatology from May 2018 is employed to represent the climatology from January to April 2013.
  - To account for prolonged daylight effects in the high-latitude regions of the northern hemisphere, observations from September, October, and November of the previous year are utilized to generate climatology. For the initial years (2012 for VNP and 2018 for VJ1) where no prior year data is available, data from the specified year (2012 for VNP and 2018 for VJ1) is used.



## Download GEOS-FP

Module **GEOS_FP.py** provide the function to download the GEOS-FP data from the GMAO.

```
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
```

To be noting, geos-fp is a forecast data generated by GMAO, for a single day, only 24 **tavg1_2d_flx_Nx** files will be generated. The latency of the **tavg1_2d_flx_Nx** can be mitigated by coordinating with GMAO office.

## FILDA2_Test_Demo

**FILDA2_Test_Demo** is a monimal test demo of FILDA2. Four directories are in the root of **FILDA2_Test_Demo**

![alt text](https://github.com/menkbard/FILDA2/blob/master/docs/FILDA_workflow.jpg)


 - FILDA2 contains all the scripts of FILDA2 package
 - Static is a mimic of the [Static data](#static-data). It worth noting that user does not need to follow the structure to put all the static data into one directory, as long as each of the static data follow the layout requirement in [Layout of the input dataset](#layout-of-the-input-dataset)
 - RUN is where the Level-1b and GOES-FP data stored. The structure shows in the image above applies to Layout_MODE=3. If COPY_STATIC=1, after excuting the FILDA2, the static data will be copied to RUN directory.
 - OUTPUT is the layout the output

## Runtime


## License
[MIT](https://choosealicense.com/licenses/mit/)

## Institution
 Atmospheric and Environmental Research Lab, the University of Iowa

## Contact

 - Principal Investigator : Prof. Jun Wang (jun-wang-1@uiowa.edu)
 - Main developer  	   : Meng Zhou (meng-zhou-1@uiowa.edu)
 - Co-developer    	   : Lorena Castro (lorena-castrogarcia@uiowa.edu)

## Website
 Visite the AER Lab website for the near-real time FILDA-2 operational
 results http://esmc.uiowa.edu:3838/fires_detection/

