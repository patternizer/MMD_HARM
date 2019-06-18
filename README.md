<img alt="FIDUCEO: MMD_HARM" align="right" src="http://www.fiduceo.eu/sites/default/files/FIDUCEO-logo.png">

# MMD_HARM

Development code for generation of a SST matchup database from AVHRR Easy FCDR + SST_CCI processing.

## Contents

* `LICENSE` - FIDCUEO license file
* `README.md` - README file
* `convert_func.py` - conversion functions for BT/radiance/counts + LUT
* `test_l1b.py` - load L1B counts and temperatures from MetOp-A orbit and calculate L using coeffs from harmonisation and convert to BT.
* `test_mmd.py` - load L1C matchup database for MetOp-A and compare with BT calculated using coeffs from harmonisation.
* `test_cci.py` - load L1C ESA SST_CCI AVHRRMTA_G v1.5 product and compare with BT calculated using coeffs from harmonisation.

## Data

Data files needed to run the tests:

* `lut_BT.npy` - standalone sensor and channel-dependent BT LUT
* `lut_radiance.npy` - standalone sensor and channel-dependent radiance LUT
* `mta_l1b.nc` - netCDF-4: Level-1B FIDUCEO AVHRRMTA_G L1B counts and temperature data for a test orbit from MetOp-A
* `mta_l1c.nc` - netCDF-4: Level-1C FIDUCEO AVHRRMTA_G L1C product v0.2Bet for a test orbit from MetOp-A
* `mta_cci.nc` - netCDF-4: Level-1C ESA SST_CCI AVHRRMTA_G L1C product v1.5 for a test orbit from MetOp-A
* `mta_mmd.nc` - netCDF-4: Level-1B FIDUCEO multi-sensor match-up dataset for a test orbit from MetOp-A 
* `FIDUCEO_Harmonisation_Data_37.nc` - netCDF-4: 3.7 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_11.nc` - netCDF-4: 11 micron channel harmonisation test data containing sensor and channel-dependent coefficients
* `FIDUCEO_Harmonisation_Data_12.nc` - netCDF-4: 12 micron channel harmonisation test data containing sensor and channel-dependent coefficients

Available on request from https://patternizer.github.io

## Contact information

* Michael Taylor (michael.taylor@reading.ac.uk)



