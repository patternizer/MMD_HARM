
<img alt="FIDUCEO: MMD_HARM" align="right" src="http://www.fiduceo.eu/sites/default/files/FIDUCEO-logo.png">

# MMD_HARM

Development code for generation of SST matchup database from AVHRR Easy FCDR + SST_CCI processing.

## Contents

* `LICENSE` - FIDCUEO license file
* `README.md` - README file
* `calc_radiance.py` - implementation of sensor and channel-dependent measurement equations (now superceded by convert_func.py)
* `convert_func.py` - conversion functions for BT/radiance/counts + LUT

## Data

Data files needed to run calc_radiance.py:

* `lut_BT.npy` - standalone sensor and channel-dependent BT LUT
* `lut_radiance.npy` - standalone sensor and channel-dependent radiance LUT
* `mta_l1b.nc` - netCDF-4: Level-1B counts and temperature data for a test orbit from MetOp-A
* `FIDUCEO_Harmonisation_Data_37.nc` - netCDF-4: 3.7 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_11.nc` - netCDF-4: 11 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_12.nc` - netCDF-4: 12 micron channel harmonisation test data containing sensor and channel-dependent coefficients

Available on request from https://github.com/patternizer

## Contact information

* Michael Taylor (michael.taylor@reading.ac.uk)



