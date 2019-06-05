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
    l1b_file = 'mta_l1b.nc'
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
    # Load: L1b orbit counts and temperatures
    #

    fcdr = xarray.open_dataset(l1b_file, decode_times=False)

    #
    # Load: harmonisation coefficients
    #

    ds = xarray.open_dataset(harm_file)
        
    if channel == 37:
        a0 = ds.parameter[(idx_ *3)].values
        a1 = ds.parameter[(idx_ *3)+1].values
        a2 = ds.parameter[(idx_ *3)+2].values
        a3 = np.nan
    else:
        a0 = ds.parameter[(idx_ *4)].values
        a1 = ds.parameter[(idx_ *4)+1].values
        a2 = ds.parameter[(idx_ *4)+2].values
        a3 = ds.parameter[(idx_ *4)+3].values

    #
    # Calculate radiance from counts and temperatures with measurement equation
    #

    if channel == 3:
        Ce = fcdr.ch3_earth_counts
        Cs = fcdr.ch3_space_counts
        Cict = fcdr.ch3_bb_counts
        Lict = fcdr.ICT_Rad_Ch3
    elif channel == 4:
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

    L = con.count2rad(Ce,Cs,Cict,Lict,Tinst,channel,a0,a1,a2,a3)

    #
    # Convert radiance to BT
    #

    BT = con.rad2bt(L,channel,lut)

    #
    # Plot orbital scan radiance and BT
    #

    bad_data = np.less_equal(BT,np.zeros(BT.shape))
    BT[bad_data] = np.nan

    fig, ax = plt.subplots(2, 1, subplot_kw={'projection': '3d'}) 
    X, Y = np.meshgrid(np.array(range(L.shape[0])), np.array(range(L.shape[1])))
    ax[0].plot_surface(X.T, Y.T, L, cmap='viridis')
    ax[1].plot_surface(X.T, Y.T, BT, cmap='viridis')
    ax[0].set_title('Radiance')
    ax[1].set_title('Brightness Temperature')
    fig.tight_layout()
    plt.savefig('L_BT.png')
    
    print('** END')
