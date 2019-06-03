import numpy as np
from netCDF4 import Dataset



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
#  count2rad()                                              #
#  REQUIRES                                                 #
#  RETURNS counts in array with same shape as L input       #
#                                                           #
#  drad_da()                                                #
#  REQUIRES                                                 #
#  RETURNS: If 3.7um channel, then two derivatives of       #
#  radiance with respect to a1 and a2 coefficents. Else     #
#  it returns three derivitives of radiance with respect    #
#  to a1, a2 and a3 coefficents                             #
#                                                           #
#############################################################





def read_in_LUT(avhrr_sat):
    LUT = {}
    all_lut_radiance_dict = np.load('lut_radiance.npy').item()
    all_lut_BT_dict = np.load('lut_BT.npy').item()
    try:
        LUT['L'] = all_lut_radiance_dict[avhrr_sat][:]
        LUT['BT'] = all_lut_BT_dict[avhrr_sat][:]
    except:
        print "Sensor for AVHRR does not exist: ", avhrr_sat
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


def count2rad(Ce,Cs,Cict,Tict,Tinst,channel,avhrr_sat):



    return counts


def drad_da(Ce,Cs,Cict,Tict,Tinst,channel,avhrr_sat):

    if channel == 37:


        return drad_da1,drad_da2
   
    else:

        return drad_da1,drad_da2,drad_da3
