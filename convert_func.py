import numpy as np
from netCDF4 import Dataset
import xarray
import matplotlib.pylab as plt


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def dbyd(a1,a2,b1,b2):
    dbyd = ((a2-a1)/(b2-b1))
    return dbyd

def load_data(file_in):
    '''
    Load netcdf into xarray
    '''
    ds = xarray.open_dataset(file_in)
    return ds

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
#  REQUIRES arrays Ce,Cs,Cict,Lict,Tinst, channel (either   #
# 37, 11 or 12), the pixel indices ni,nj, and the           #
# harmonisation coefficients a0,a1,a2,a3                    #
#  RETURNS radiance in array with same shape as L input     #
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

def drad_da(Ce,Cs,Cict,Tict,Tinst,channel,avhrr_sat):
    if channel == 37:
        return drad_da1,drad_da2
    else:
        return drad_da1,drad_da2,drad_da3

def count2rad(Ce,Cs,Cict,Lict,Tinst,channel,avhrr_sat,ni,nj,a0,a1,a2,a3):
    L = np.empty(shape=(len(nj),len(ni)))
    e0 = 0.985140
    if channel == 37:
#        L = (a0 + ((Lict * (e0 + a1)) / (Cict - Cs)) * (Ce - Cs) + a2 * Tinst) * 100.0
        L = a0 + ((Lict * (e0 + a1)) / (Cict - Cs)) * (Ce - Cs) + a2 * Tinst
    else:
        L = a0 + ((Lict * (e0 + a1)) / (Cict - Cs) + a2 * (Ce - Cict)) * (Ce - Cs) + a3 * Tinst
    return L

# =======================================
# RUN TEST CASE
# =======================================

if __name__ == "__main__":

    avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']
    T_mean = np.array([ 290.327113, 288.637636, 294.758564, 292.672201, 288.106630, 288.219774, 287.754638, 286.125823, np.nan])
    T_sdev = np.array([ 2.120666, 1.053762, 2.804361, 3.805704, 1.607656, 0.607697, 0.117681, 0.049088, np.nan])

    #
    # Select channel
    #

    channel = 37

    #
    # Read in lut for MTA
    #

    idx = 7
    lut = read_in_LUT(avhrr_sat[idx]) # MTA

    #
    # Load L1b orbit file from MTA
    #

    fcdr = xarray.open_dataset("mta_l1b.nc", decode_times=False)

    #
    # Load harmonisation coefficients
    #

    file_in = "FIDUCEO_Harmonisation_Data_" + str(channel) + ".nc"
    ds = load_data(file_in)

    # avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']
    #                0       1      2      3      4      5      6     7       8
    #                MTA,  N19,   N18,   N17,    N16,   N15,  N14,    N12,   N11
    #                0       1      2      3      4      5      6     7       8

    isensor = 7 - idx
        
    if channel == 37:

        a0 = ds.parameter[(isensor*3)].values
        a1 = ds.parameter[(isensor*3)+1].values
        a2 = ds.parameter[(isensor*3)+2].values
        a3 = np.nan

    else:

        a0 = ds.parameter[(isensor*4)].values
        a1 = ds.parameter[(isensor*4)+1].values
        a2 = ds.parameter[(isensor*4)+2].values
        a3 = ds.parameter[(isensor*4)+3].values

    #
    # Calculation radiance from counts using measurement equations
    #

    ni = fcdr.nx
    nj = fcdr.ny
    if channel == 37:

        Ce = fcdr.ch3_earth_counts
        Cs = fcdr.ch3_space_counts
        Cict = fcdr.ch3_bb_counts
        Lict = fcdr.ICT_Rad_Ch3

    elif channel == 11:

        Ce = fcdr.ch4_earth_counts
        Cs = fcdr.ch4_space_counts
        Cict = fcdr.ch4_bb_counts
        Lict = fcdr.ICT_Rad_Ch4

    else:

        Ce = fcdr.ch5_earth_counts
        Cs = fcdr.ch5_space_counts
        Cict = fcdr.ch5_bb_counts
        Lict = fcdr.ICT_Rad_Ch5

    Tinst = (fcdr.prt - T_mean[idx]) / T_sdev[idx]
    L = count2rad(Ce,Cs,Cict,Lict,Tinst,channel,avhrr_sat,ni,nj,a0,a1,a2,a3)

    #
    # Convert radiance to BT
    #

    BT = rad2bt(L,channel,lut)

    #
    # Plot along track radiance and BT
    #

    mid_scan = int(len(fcdr.nx)/2)
    bad_data = BT[:,mid_scan] < 250
    BT[bad_data,mid_scan] = np.nan

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(L[:,mid_scan], 'r')
    ax[1].plot(BT[:,mid_scan], 'k')
    ax[0].set_title('Radiance')
    ax[1].set_title('Brightness Temperature')
    fig.tight_layout()
    plt.savefig('L_BT.png')
    
    #
    # TEST: read in L & BT and then do 2-way pass: rad2bt, bt2rad 
    #
    
    L = lut['L']
    BT = lut['BT']
    BT2 = rad2bt(L, channel, lut)
    L2 = bt2rad(BT2, channel, lut)
    if channel == 37:
        ch_index = 3
    elif channel == 11:
        ch_index = 4
    else:
        ch_index = 5        

    fig, ax = plt.subplots(2,1)
    ax[0].plot(L2[:,ch_index],'.')
    ax[0].plot(L[:,ch_index],'-')
    ax[1].plot(BT2[:,ch_index],'.')
    ax[1].plot(BT[:,ch_index],'-')
    ax[0].set_title('Radiance')
    ax[1].set_title('Brightness Temperature')
    fig.tight_layout()
    pltstr = 'test_ch_' + str(ch_index) + '.png'
    plt.savefig(pltstr)

    print('** END')

