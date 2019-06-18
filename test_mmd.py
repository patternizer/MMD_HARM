# =======================================
# RUN TEST CASE
# =======================================
import numpy as np
from netCDF4 import Dataset
import xarray
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import axes3d
import convert_func as con

if __name__ == "__main__":

    ################################################################
    # RUN PARAMETERS: 
    ################################################################

    ch = 37
    idx = 7      # MTA (see avhrr_sat)
    mmd_file = 'mta_mmd.nc'
    harm_file = 'FIDUCEO_Harmonisation_Data_' + str(ch) + '.nc'
    if ch == 37:
        channel = 3
    elif ch == 11:
        channel = 4
    else:
        channel = 5

    ################################################################

    #
    # senesor list, Tinst normalisation parameters: T_mean and T_sdev
    #

    avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']

    # NB: Harmonisation coefficients provided by Ralf Quast are for MTA --> N11:
    # RQ:          MTA,   N19,   N18,   N17,   N16,   N15,  N14,    N12,   N11
    # index         0       1      2      3      4      5      6     7       8
    # --> new index map

    idx_ = 7 - idx

    T_mean = np.array([ 290.327113, 288.637636, 294.758564, 292.672201, 288.106630, 288.219774, 287.754638, 286.125823, np.nan])
    T_sdev = np.array([ 2.120666, 1.053762, 2.804361, 3.805704, 1.607656, 0.607697, 0.117681, 0.049088, np.nan])

    #
    # Load: lut
    #

    lut = con.read_in_LUT(avhrr_sat[idx])

    #
    # Load: harmonisation coefficients
    #

    ds = xarray.open_dataset(harm_file)
        
    if channel == 37:
        a1 = ds.parameter[(idx_ *3)].values
        a2 = ds.parameter[(idx_ *3)+1].values
        a3 = ds.parameter[(idx_ *3)+2].values
        a4 = 0.0
        a5 = 0.0
    else:
        a1 = ds.parameter[(idx_ *4)].values
        a2 = ds.parameter[(idx_ *4)+1].values
        a3 = ds.parameter[(idx_ *4)+2].values
        a4 = ds.parameter[(idx_ *4)+3].values
        a5 = 0.0

    #
    # Load: L1b orbit counts and temperatures
    #

    mmd = xarray.open_dataset(mmd_file, decode_times=False)

    if channel == 3:
        Ce = mmd['avhrr-ma_ch3b_earth_counts']
        Cs = mmd['avhrr-ma_ch3b_space_counts']
        Cict = mmd['avhrr-ma_ch3b_bbody_counts']
        # Lict = fcdr['ICT_Rad_Ch3']

    elif channel == 4:
        Ce = mmd['avhrr-ma_ch4_earth_counts']
        Cs = mmd['avhrr-ma_ch4_space_counts']
        Cict = mmd['avhrr-ma_ch4_bbody_counts']
        # Lict = fcdr['ICT_Rad_Ch4']

    else:
        Ce = mmd['avhrr-ma_ch5_earth_counts']
        Cs = mmd['avhrr-ma_ch5_space_counts']
        Cict = mmd['avhrr-ma_ch5_bbody_counts']
        # Lict = fcdr['ICT_Rad_Ch5']

#
#   convert Tict --> Lict
#

    Tict = mmd['avhrr-ma_ict_temp']
    Lict = con.bt2rad(Tict,channel,lut)

    Tinst = (mmd['avhrr-ma_orbital_temperature'] - T_mean[idx]) / T_sdev[idx]
    WV = 0.0 * Tinst

    #
    # Calculate radiance from counts and temperatures with measurement equation
    #

    L = con.count2rad(Ce,Cs,Cict,Lict,Tinst,WV,channel,a1,a2,a3,a4,a5)

    #
    # Convert radiance to BT
    #

    BT = con.rad2bt(L,channel,lut)

    #
    # Plot BT versus BT in MMD at centre pixel (3,3)
    #

    if channel == 3: BT_MMD = mmd['avhrr-ma_ch3b'][:,3,3]
    elif channel == 4: BT_MMD = mmd['avhrr-ma_ch4'][:,3,3]
    else: BT_MMD = mmd['avhrr-ma_ch5'][:,3,3]
    
    BT = BT[:,3,3]
    gd = BT > 0

    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(BT[gd],'.',markersize=0.5, label='BT (from ME)')
    ax[0,0].plot(BT_MMD[gd],'.',markersize=0.5, label='BT (from MMD)')
    ax[0,0].legend(markerscale=20, scatterpoints=5, fontsize=8)
    ax[0,0].set_ylim(200,320)
    ax[0,0].set_xlabel('MM')
    ax[0,0].set_ylabel('BT')
    ax[1,0].plot(BT_MMD[gd], BT_MMD[gd]-BT[gd],'.',markersize=0.5)
    ax[1,0].set_xlim(200,320)
    ax[1,0].set_ylim(-1,1)
    ax[1,0].set_xlabel('BT_MMD')
    ax[1,0].set_ylabel('BT(MMD)-BT(ME)')
#    ax[0,1].hist((BT[gd]-float(BT[gd].mean())),bins=100, label='d_BT (from ME)')
#    ax[0,1].hist((BT_MMD[gd]-float(BT_MMD[gd].mean())),bins=100, label='d_BT (from MMD)')
#    ax[0,1].legend(markerscale=20, scatterpoints=5, fontsize=8)
    ax[0,1].hist((BT_MMD[gd]-BT[gd]),bins=100)
    ax[0,1].set_xlabel('BT(MMD)-BT(ME)')
    ax[1,1].set_ylabel('count (nbins=100)')
    ax[1,1].plot(BT_MMD[gd], BT_MMD[gd], '.',markersize=0.5)
    ax[1,1].set_xlim(200,320)
    ax[1,1].set_ylim(200,320)
    ax[1,1].set_xlabel('BT_MMD')
    ax[1,1].set_ylabel('BT(ME)')
    plt.tight_layout()
    plt.savefig('BT_v_BT_MMD.png')
    
    print('** END')
