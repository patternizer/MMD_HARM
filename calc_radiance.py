#!/usr/bin/env python

# call as: python calc_radiance.py

# =======================================
# Version 0.1
# 31 May, 2019
# michael.taylor AT reading DOT ac DOT uk
# =======================================

import os  
import os.path  
from os import fsync, remove
import glob  
import optparse 
from  optparse import OptionParser  
import sys   
import numpy as np
import numpy.ma as ma  
from numpy import array_equal, savetxt, loadtxt, frombuffer, save as np_save, load as np_load, savez_compressed, array
import xarray
import pandas as pd 
from pandas import Series, DataFrame, Panel 
import matplotlib.pyplot as plt

# =======================================
# METHODS
# =======================================    

def load_data(file_in):
    '''
    Load netcdf into xarray
    '''
    ds = xarray.open_dataset(file_in)
    return ds

def calc_radiance(ds, isensor, nch, fcdr):
    '''
    Calculate radiance using channel-dependent measurement equations
    '''

    ni = fcdr.nx
    nj = fcdr.ny
    C_ict_37 = fcdr.ch3_bb_counts  
    C_s_37 = fcdr.ch3_space_counts
    C_e_37 = fcdr.ch3_earth_counts
    L_ict_37 = fcdr.ICT_Rad_Ch3    
    C_ict_11 = fcdr.ch4_bb_counts  
    C_s_11 = fcdr.ch4_space_counts 
    C_e_11 = fcdr.ch4_earth_counts 
    L_ict_11 = fcdr.ICT_Rad_Ch4    
    C_ict_12 = fcdr.ch5_bb_counts  
    C_s_12 = fcdr.ch5_space_counts 
    C_e_12 = fcdr.ch5_earth_counts 
    L_ict_12 = fcdr.ICT_Rad_Ch5   
    T_inst = fcdr.prt              
    
    e_ict = 0.985140

    # | Sensor | T_min (K) | T_max (K) | T_mean (K) | T_sdev (K) | Meas eqn (11, 12 / 3.7 µm) |
    # |--------|-----------|-----------|------------|------------|----------------------------| 
    # | m02    | 285.9     | 286.4     | 286.125823 | 0.049088   | 102 / 106                  |
    # | n19    | 286.9     | 288.3     | 287.754638 | 0.117681   | 102 / 106                  |
    # | n18    | 286.6     | 290.2     | 288.219774 | 0.607697   | 102 / 106                  |
    # | n17    | 286.1     | 298.2     | 288.106630 | 1.607656   | 102 / 106                  |
    # | n16    | 287.2     | 302.0     | 292.672201 | 3.805704   | 102 / 106                  |
    # | n15    | 285.1     | 300.6     | 294.758564 | 2.804361   | 102 / 106                  |
    # | n14    | 286.8     | 296.4     | 288.637636 | 1.053762   | 102 / 106                  |
    # | n12    | 287.2     | 302.8     | 290.327113 | 2.120666   | 102 / 106                  |
    # | n11    | 286.1     | 299.9     | 290.402168 | 3.694937   | 102 / 106                  |
  
    T_ave = np.array([ 286.125823, 287.754638, 288.219774, 288.106630, 292.672201, 294.758564, 288.637636, 290.327113, 290.402168 ])
    T_std = np.array([ 0.049088, 0.117681, 0.607697, 1.607656, 3.805704, 2.804361, 1.053762, 2.120666, 3.694937])
    T_mean = T_ave[isensor]
    T_sdev = T_std[isensor]

    L = np.empty(shape=(len(nj),len(ni)))    
    if nch == 37:

        a0 = ds.parameter[(isensor*3)]
        a1 = ds.parameter[(isensor*3)+1]
        a2 = ds.parameter[(isensor*3)+2]

        # Measurement equation 106:
        L = a0 + ((L_ict_37 * (e_ict + a1)) / (C_ict_37 - C_s_37)) * (C_e_37 - C_s_37) + a2 * (T_inst - T_mean) / T_sdev

        # NB: for 3.7 micron channel, scale by factor of 100 to get correct radiance
        L = L * 100.0

    elif nch == 11:

        a0 = ds.parameter[(isensor*4)]
        a1 = ds.parameter[(isensor*4)+1]
        a2 = ds.parameter[(isensor*4)+2]
        a3 = ds.parameter[(isensor*4)+3]

        # Measurement equation 102:
        L = a0 + ((L_ict_11 * (e_ict + a1)) / (C_ict_11 - C_s_11) + a2 * (C_e_11 - C_ict_11)) * (C_e_11 - C_s_11) + a3 * (T_inst - T_mean) / T_sdev

    else:

        a0 = ds.parameter[(isensor*4)]
        a1 = ds.parameter[(isensor*4)+1]
        a2 = ds.parameter[(isensor*4)+2]
        a3 = ds.parameter[(isensor*4)+3]

        # Measurement equation 102:
        L = a0 + ((L_ict_12 * (e_ict + a1)) / (C_ict_12 - C_s_12) + a2 * (C_e_12 - C_ict_12)) * (C_e_12 - C_s_12) + a3 * (T_inst - T_mean) / T_sdev

    return L

def radiance2bt(L, nch, lut):
    '''
    Use sensor-specific look-up tables to convert radiance (L) to brightness temperature (BT)
    '''

    lut_L = lut.lookup_table_radiance
    lut_BT = lut.lookup_table_BT

    BT = np.empty(shape=(L.shape[0],L.shape[1]))

    if nch == 37: channel = 3
    elif nch == 11: channel = 4
    else: channel = 5

    BT = np.interp(L, lut_L[:,channel], lut_BT[:,channel])

    return BT

def bt2radiance(BT, nch, lut):
    '''
    Use sensor-specific look-up tables to convert brightness temperature (BT) to radiance (L)
    '''

    lut_L = lut.lookup_table_radiance
    lut_BT = lut.lookup_table_BT

    L = np.empty(shape=(BT.shape[0],BT.shape[1]))

    if nch == 37: channel = 3
    elif nch == 11: channel = 4
    else: channel = 5

    L = np.interp(BT, lut_BT[:,channel], lut_L[:,channel])

    return L
        
# =======================================    
# MAIN BLOCK
# =======================================    
    
if __name__ == "__main__":

    #--------------------------------------------------
    # parser = OptionParser("usage: %prog nch isensor")
    # (options, args) = parser.parse_args()

    # nch = args[0]
    # isensor = args[1]

    nch = 37
    # nch = 11
    # nch = 12

    sensor = ['METOPA','NOAA19','NOAA18','NOAA17','NOAA16','NOAA15','NOAA14','NOAA12','NOAA11']
    isensor = 0

    #--------------------------------------------------

    #
    # Load harmonisation file
    #

    file_in = "FIDUCEO_Harmonisation_Data_" + str(nch) + ".nc"
    ds = load_data(file_in)

    #
    # Load L1B orbit FCDR data (fcdr) and radiance <--> BT look-up table (lut)
    #

    if isensor == 0:

        # MetOp-A:

        fcdr = xarray.open_dataset("mta_l1b.nc", decode_times=False) 
        lut = xarray.open_dataset("mta_l1c.nc")
        L = calc_radiance(ds, isensor, nch, fcdr)
        BT = radiance2bt(L, nch, lut)
        L2 = bt2radiance(BT, nch, lut)

        mid_scan = int(len(fcdr.nx)/2)

        fig, ax = plt.subplots(2, 1)
        ax[0].plot(L[:,mid_scan], 'r')
        ax[1].plot(BT[:,mid_scan], 'k')
        ax[0].set_title('Radiance')
        ax[1].set_title('Brightness Temperature')
        fig.tight_layout()
        plt.savefig('L_BT.png')

    else:
        
        print('*** No sensor chosen')

    print('*** End of program')


