# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Deprecated

## [v1.0 2024-10-08]

### Added
- Add Option to determine the number of core in MP automatically
- Add VIIRS M10, M08 sensor data for the estimation
- Add a science value **P lower bound** to specify the minimal values for the fire fraction
- Add output of FP_Power_R_AC. Radiance based FRP estimation with atmospheric correction

### Changed
- Change the license of MCBEF package to Apache 2.0
- Redefine the prior for the uni-phasic model
- Redefine the selection of fire model for the estimation
	- McBEF now will first attempt to use bi-phasic model for parameter estimates for all the input regardless of the **Gas flaring** or **Static source** flag. If fail, it will then try the uni-phasic model for the estimation


### Fixed

### Removed
- Old MBFPE package

### Deprecated
