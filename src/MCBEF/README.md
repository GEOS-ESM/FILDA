


# McBEF
Monte Carlo Biphasic Estimation of Fire Properties (McBEF)

The McBEF is an algorithm implementation to utilizes multiple fire-sensitive 
channels from visible to thermal infrared, available on advanced sensors 
like the Visible Infrared Imaging Radiometer Suite (VIIRS) and Geostationary
 Operational Environmental Satellite (GOES) for sub-pixel fire flaming and smoldering
  temperature estimation.
  
##  - **Recent update** - 

1. Add VIIRS M10, M08 sensor data for the estimation
2. Add a science value **P lower bound** to specify the minimal values for the fire fraction
3. Add output of FP_Power_R_AC. Radiance based FRP estimation with atmospheric correction
4. Redefine the prior for the uni-phasic model
5. Redefine the selection of fire model for the estimation
	- McBEF now will first attempt to use bi-phasic model for parameter estimates for all the input regardless of the **Gas flaring** or **Static source** flag. If fail, it will then try the uni-phasic model for the estimation.

	

## Table of content


- [Package Dependence](#package-dependence)

- [How to use McBEF](#how-to-use-McBEF)
	- [Forward processing](#Forward-processing)
	- [Reprocessing](#Reprocessing)

## Package Dependence

| #  |    libraries     | version   | Description| 
| ---|       :---:      |   :---:   | ------------- |
| 1  | os               | Default   | Miscellaneous operating system interfaces |
| 2  | sys              | Default   | System-specific parameters and functions  |
| 3  | copy             | Default   | Shallow and deep copy operations  |
| 4  | datetime         | Default   | Basic date and time types  |
| 5  | multiprocessing  | Default   | Process-based parallelism  |
| 6  | numpy            | 1.21.5    | https://numpy.org/  |
| 7  | scipy            | 1.7.3     | https://scipy.org/  |
| 8  | pandas           | 1.3.5     | https://pandas.pydata.org/	  |
| 9  | netCDF4          | 1.5.7     | https://unidata.github.io/netcdf4-python/  |
| 10 | h5py             | 3.7.0     | https://docs.h5py.org/en/stable/build.html#source-installation |
| 11 | pymc3            | 3.11.2    | https://www.pymc.io/projects/examples/en/2021.11.0/getting_started.html |
| 11 | theano-pymc      | 1.1.2     | https://pypi.org/project/Theano-PyMC/ |
| 11 | arviz            | 0.11.2    | https://www.arviz.org/en/latest/ |


**NOTE**:

The McBEF was originally developed using the PyMC3 framework, which is the 3.0 version of the PyMC project, utilizing Theano as its tensor backend.

The PyMC development team later rebranded the project as PyMC and replaced Theano with Aesara as the tensor backend. More details can be found in the official documentation:

[PyMC v4 Announcement](https://www.pymc.io/blog/v4_announcement.html#v4_announcement)

As of now, the latest PyMC release is version 5.16.2, which uses PyTensor as the backend. Relevant documentation can be accessed here:

[PyMC API Documentation](https://www.pymc.io/projects/docs/en/stable/api.html)

[PyTensor Overview](https://www.pymc.io/projects/docs/en/stable/learn/core_notebooks/pymc_pytensor.html)

After testing nearly all versions of PyMC with their respective tensor libraries across several case studies, I did not observe significant differences in estimation results (though not in speed). Given the somewhat confusing documentation on the transitions between Theano, Aesara, and PyTensor, I have opted to stick with the vanilla version of PyMC for now until a more stable and well-documented version is released.

```FILDA.yml``` contains the list of the necessary libraries to run FILDA-2. We recommend to create an environment and
install all of these libraries at the time. The above can be performed by using ```Anaconda``` or ```Pip``` as follows.

### Anaconda
If you work under Anaconda (preferable, a library manage tool that help to check the interdenpendency of libraries), you can create an anaconda 
environment to install all the libraries by using the following anaconda command:
```
conda env create -f MCBEF.yml
```

To run McBEF in the McBEF python environment, you need to activate the environment using the following command:

```
source activate MCBEF
```
### Pip
Alternatively, the list of McBEF requirements can be installed into a standard python virtual environment and then 
installing the libraries by using the ```pip``` packages insalled. The general steps for this approach are:
1. Create a new environment
```
    python -m venv env
```
2. Activate the virtual environment
```
    source env/bin/activate
```
3. Install the libraries specified in the MCBEF.yml file
```
    pip install -r MCBEF.yml
```

In both cases, once the environment is created and the requirement are installed, you should be able to implement the McBEF library.


## How to use McBEF:

### Forward processing

The McBEF process the VNP47IMG Level2 prdouct at overpass level, as input parameters they are required:
```
- jdn, julian day in the format "AYYYYDDD", e.g. "A2022006"
- overpass_beg, start time in the format "HHMM", e.g. "0018"
- overpass_end, stop time in the format "HHMM", e.g. "0024"
- sat, satellite, options: VNP | VJ1
```
**NOTE**:
Overpass are generated within the half-open interval [overpass_beg, overpass_end) (in other words, the interval including start time but excluding stop time).

For case of using a *bash script* to control input parameters for running the package, use script in the following manners:
```
python main.py A2019213 0000 2400 VNP
```

### Reprocessing
```
python main_reprocessing.py A2019213 0000 2400 VNP
```

## Namelist
The namelist is the configuration file to control the McBEF package. There are three blocks in the current namelist, **CONTROL MENU**, **I/O MENU**, **SCIENCE VALUE**. A operational user will only need to mofidy the first and second section. 

### CONTROL MENU
```
# BPFPE V1.0 namelist file Aug. 2024, Meng Zhou, Arlindo da Silva
------------------------+-----------------------------------------------
%%%   CONTROL MENU   %%%:
Background Band         : M14 M15 I05 M16
Fire Band               : DNB M11 M13 M14 M15 M16
Number of Draw          : 1000
Number of Tune          : 2000
Number of Chain         : 2
Number of Core          : 70
Prior distribution      : U
Emissivity NRT          : T
Gradient Sampling?      : F
Turn on debug printout? : F
```
 - **Background Band**: Band observations used to estimate the background temperature and the water vapor scaling facor, following the naming convention of the VIIRS instruments
 - **Fire Band**: Band observations used to estimate the fire parameter
 - **Number of Draw**: Number of samples that will be collected from the posterior distribution for each chain
 - **Number of Tune**: Number of tuning steps that the sampler will execute before beginning to collect the actual samples.
 - **Number of Chain**: Number of independent Markov Chains that samples will be drawn from.
 - **Number of Core**: Number of cores that allows for multi-processing
	 - if **-1** is provided, McBEF will get the number of cores within the node (N) and use N-4 as the available cores for multiprocessing
 - **Prior distribution**: Options for the prior distribution, U for uniform distribution, G for Gaussian distribution 
 - **Emissivity NRT**: If **True**, McBEF will try to find the near-real time VIIRS/NPP Land Surface Temperature (VNP21A2).
 - **Gradient Sampling**: If **True**, McBEF will use the No U-Turn Sampler ([NUTS](https://www.pymc.io/projects/docs/en/stable/api/generated/pymc.step_methods.hmc.NUTS.html)).
 - **Turn on debug printout?**: If **True**, McBEF will printout intermediate results. 
### I/O MENU
```
# BPFPE V1.0 namelist file Aug. 2024, Meng Zhou, Arlindo da Silva
------------------------+-----------------------------------------------
%%%      I/O MENU    %%%:
Root data directory     : ./DATA/
Run directory           : ./
Output directory        : ./
Output version          : v1_0_0
```
Some specific layout of input dataset is required, please refer to [Layout of the input dataset](#layout-of-the-input-dataset) for the preparation of the static data. 

## Input data

BPFPE relays on several inputs for fire paramter estimation, they are different in the time latency.

**1. VIIRS FILDA2 Fire Modified Combustion Efficiency products (VNP47MOD and VJ147MOD).** The FILDA-2 data can be retrieved from LAADS data archive https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP47MOD/

**2. VIIRS/NPP Land Surface Temperature/Emissivity 8-Day L3 Global 1 km (VNP21A2)** This product is updated by every 8 Days. This data can be retrieved from LAADS data archive https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP21A2/

### Climatology and Static data
**3. VIIRS/NPP Land Surface Temperature/Emissivity L3 Global 1 km  climatology**

**4. Static thermal anomaly**
The Static thermal anomaly dataset is a rasterized dataset. The format of the dataset is netCDF4 (307K).

**=====================================**
**IMPORTANT NOTE for dataset**
**=====================================**
The BPFPE package adopts a flexible approach in selecting the required datasets, such as VNP21A2. These datasets are chosen based on availability but are not mandatory for the process. In essence, BPFPE evaluates the necessary files for processing a VIIRS granule and attempts to locate them within the provided directory specified in the namelist. If BPFPE successfully finds a file, it will read and utilize it accordingly. However, if a file is not found, it will raise a warning instead of breaking the code or triggering an error. This is particularly relevant for ocean surface data, as there may be cases where no dataset is available from the original provider.

## Layout of the input dataset
Below is layout of the input DATA directory, user need to follow this structure to allow BPFPE find those necessary input data.

- DATA
  - STATIC
    - static_thermal_anomaly.2019.nc
  - VNP21A2
    - 2013
    - ...
    - 2023
        - 001 
        - ...
        - 361
    - CLT: Climatology data generated by three year average of the VNP21A2
        - 001 
        - ...
        - 361: Data under each day of year are arranged by Sinusodial tiles
            - VNP21A2.A2019009.h00v08.clt.nc
            - ...
  - VNP47MOD
    - 2013
    - ...
    - 2019
        - 001
        - ...
        - 365: Level-2 data under each day of year are arranged by severy 6 minutes
            - VNP47MOD.A2019365.0000.002.2024265020032.nc
            - ...
 - **REPROCESSING**
	 - DATA: Preprocessed daily data
		 - Intermediate.VNP.A2019183.0000.2400.nc 
	 - FIRE_MASK

## Runtime

## License
[Apache License Version 2.0](https://choosealicense.com/licenses/apache-2.0/)

## Institution

Global Modeling and Assimilation Office (GMAO), the National Aeronautics and Space Administration (NASA)

Goddard Earth Sciences Technology and Research (GESTAR) II, University of Maryland, Baltimore County

## Contact
 * Main developer : Meng Zhou (MZ)
 * Supervisor     : Arlindo da Silva (ADS)