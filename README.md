
<img alt="FIDUCEO: MMD_HARM" align="right" src="http://www.fiduceo.eu/sites/default/files/FIDUCEO-logo.png">

# MMD_HARM

Development code for generation of SST matchup database from AVHRR Easy FCDR + SST_CCI processing.

## Contents

* `LICENSE` - FIDCUEO license file
* `README.md` - README file
* `calc_radiance.py` - read in harmonisation coefficients + L1B counts data + L1C look-up tables and use measurement equations to generate radiances
* `convert_func.py` - conversion functions for BT/radiance/counts + LUT
* `lut_BT.npy` - standalone sensor and channel-dependent BT LUT
* `lut_radiance.npy` - standalone sensor and channel-dependent radiance LUT

## Data

Data files needed to run calc_radiance.py:

* `mta_l1b.nc` - netCDF-4: Level-1B counts and temperature data for a test orbit from MetOp-A
* `mta_l1c.nc` - netCDF-4: Level-1C Easy FCDR test orbit containing radiance and BT look-up tables from MetOp-A
* `FIDUCEO_Harmonisation_Data_37.nc` - netCDF-4: 3.7 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_11.nc` - netCDF-4: 11 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_12.nc` - netCDF-4: 12 micron channel harmonisation test data containing sensor and channel-dependent coefficients

Available on request from https://github.com/patternizer

## Contact information

* Michael Taylor (michael.taylor@reading.ac.uk)



