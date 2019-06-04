import numpy as np
from netCDF4 import Dataset
import xarray
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import axes3d

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def dbyd(a1,a2,b1,b2):
    dbyd = ((a2-a1)/(b2-b1))
    return dbyd

#############################################################
#                                                           #
#  FUNCTIONS FOR CONVERTING BETWEEN BT, COUNTS AND RADIANCE #
#                                                           #
#  read_in_LUT()                                            #
#  REQUIRES avhrr_sat to be in the format:                  #
#  N12,N14,N15,N16,N17,N18,N19,MTA,MTB                      #
#  RETURNS LUT dictionary with groups "L" and "BT"          #
#                                                           #
#  rad2bt()                                                 #
#  REQUIRES array of radiance, channel (either 37, 11 or 12)#
#  and lut outputfrom read_in_LUT (sensor specific)         #
#  RETURNS BT in array of same initial shape of L           #
#                                                           #
#  bt2rad()                                                 #
#  REQUIRES array of BT, channel (either 37, 11 or 12)      #
#  and lut output from read_in_LUT (sensor specific)        #
#  RETURNS L in array of same initial shape of BT           #
#                                                           #
#  dbt_drad()                                               #
#  REQUIRES array of radiance, channel (either 37, 11 or 12)#
#  and lut output from read_in_LUT (sensor specific)        #
#  RETURNS dbt_drad array in the same shape as L input      #
#                                                           #
#  drad_da()                                                #
#  REQUIRES                                                 #
#  RETURNS: If 3.7um channel, then two derivatives of       #
#  radiance with respect to a1 and a2 coefficents. Else     #
#  it returns three derivitives of radiance with respect    #
#  to a1, a2 and a3 coefficents                             #
#                                                           #
#  count2rad()                                              #
#  REQUIRES arrays Ce,Cs,Cict,Lict,Tinst, channel (either   #
#  37, 11 or 12), and harmon coefficients a0,a1,a2,a3       #
#  RETURNS radiance in array with same shape as L input     #
#                                                           #
#############################################################

def read_in_LUT(avhrr_sat):
    LUT = {}
    all_lut_radiance_dict = np.load('lut_radiance.npy', encoding='bytes').item()
    all_lut_BT_dict = np.load('lut_BT.npy', encoding='bytes').item()
    try:
        LUT['L'] = all_lut_radiance_dict[avhrr_sat][:]
        LUT['BT'] = all_lut_BT_dict[avhrr_sat][:]
    except:
#        print "Sensor for AVHRR does not exist: ", avhrr_sat
        print("Sensor for AVHRR does not exist: ", avhrr_sat)
    return LUT

def rad2bt(L,channel,lut):
    if channel == 37:
        ch_index = 3
    elif channel == 11:
        ch_index = 4
    else:
        ch_index = 5
    BT = np.interp(L,lut['L'][:,ch_index],lut['BT'][:,ch_index],left=-999.9,right=-999.9)
    return BT

def bt2rad(bt,channel,lut):
    if channel == 37:
        ch_index = 3
    elif channel == 11:
        ch_index = 4
    else:
        ch_index = 5
    L = np.interp(BT,lut['BT'][:,ch_index],lut['L'][:,ch_index],left=-999.9,right=-999.9)
    return L

def dbt_drad(L,channel,lut):
    #determine channel index needed
    if channel == 37:
        ch_index = 3
    elif channel == 11:
        ch_index = 4
    else:
        ch_index = 5
    dbtdrad = np.zeros_like(L)
    # Over array L, determine the two values in the LUT either side
    # Then, find dbt by drad and add to new array
    for i in xrange(0,len(L)):
        for j in xrange(0,len(L[0])):
            element = L[i,j]
            idx = find_nearest(lut['L'][:,ch_index],element)
            if lut['L'][idx,ch_index] > element:
                dbtdrad[i,j] = (dbyd(lut['BT'][idx-1,ch_index],lut['BT'][idx,ch_index],lut['L'][idx-1,ch_index],lut['L'][idx,ch_index]))
            else:
                dbtdrad[i,j] = (dbyd(lut['BT'][idx,ch_index],lut['BT'][idx+1,ch_index],lut['L'][idx,ch_index],lut['L'][idx+1,ch_index]))
    return dbtdrad



def drad_da(L_ict,C_e,C_s,C_ict,T_ict,T_inst,T_mean,T_sdev,channel,avhrr_sat):
    # equation 106 for channel 3.7
    if channel == 37:
        drad_da1 = L_ict/(C_ict - C_s)
        drad_da2 = (T_inst - T_mean) / T_sdev
        return drad_da1,drad_da2
    # equation 102 for 11/12 channel
    else:
        drad_da1 = L_ict/(C_ict - C_s)
        drad_da2 = (C_e - C_ict) * (C_e - C_s)
        drad_da3 = (T_inst - T_mean) / T_sdev
        return drad_da1,drad_da2,drad_da3
    
    

def count2rad(Ce,Cs,Cict,Lict,Tinst,channel,a0,a1,a2,a3):
    L = np.empty(shape=(Ce.shape[0],Ce.shape[1]))
    if channel == 37:
        L = a0 + ((Lict * (0.985140 + a1)) / (Cict - Cs)) * (Ce - Cs) + a2 * Tinst
    else:
        L = a0 + ((Lict * (0.985140 + a1)) / (Cict - Cs) + a2 * (Ce - Cict)) * (Ce - Cs) + a3 * Tinst
    return L


