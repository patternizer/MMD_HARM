#!/usr/bin/env python

# Compare L1C MMD and ME brightness temperatures
# ======================================================
# Version 0.7
# 21 June, 2019
# https://patternizer.github.io
# michael.taylor AT reading DOT ac DOT uk
# ======================================================


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

    flag_new = False # NEW harmonisation structure (run >= '3.0-4d111a1')
    ch = 37
    idx = 7         # MTA (see avhrr_sat)
    if idx == 7:
        noT = True
    else:
        noT = False
    mmd_file = 'mta_mmd.nc'
    har_file = 'FIDUCEO_Harmonisation_Data_' + str(ch) + '.nc'
    if ch == 37:
        channel = 3
    elif ch == 11:
        channel = 4
    else:
        channel = 5

    ################################################################

    #
    # senesor list and indexing
    #

    # NB: LUT sensor sequence

    avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']

    # NB: HAR sensor sequence

    if flag_new:
        # RQ: NEW   AATSR,   MTA,   N19,   N18,   N17,   N16,   N15,   N14,   N12,   N11
        # index         0      1      2      3      4      5      6      7      8      9
        # --> new index map
        idx_ = 7 - idx + 1 
    else:
        # RQ: OLD     MTA,   N19,   N18,   N17,   N16,   N15,   N14,   N12,   N11
        # index         0      1      2      3      4      5      6      7      8 
        # --> new index map
        idx_ = 7 - idx

    #
    # Load: lut
    #

    lut = con.read_in_LUT(avhrr_sat[idx])

    #
    # Load: mmd orbit counts and temperatures
    #

    mmd = xarray.open_dataset(mmd_file, decode_times=False)

    if channel == 3: 
        BT_MMD = mmd['avhrr-ma_ch3b']    # (55604, 7, 7)
    elif channel == 4: 
        BT_MMD = mmd['avhrr-ma_ch4']
    else: 
        BT_MMD = mmd['avhrr-ma_ch5']
     
    #
    # Load: mmd counts and temperatures
    #

    if channel == 3:
        Ce = mmd['avhrr-ma_ch3b_earth_counts']
        Cs = mmd['avhrr-ma_ch3b_space_counts']
        Cict = mmd['avhrr-ma_ch3b_bbody_counts']

    elif channel == 4:
        Ce = mmd['avhrr-ma_ch4_earth_counts']
        Cs = mmd['avhrr-ma_ch4_space_counts']
        Cict = mmd['avhrr-ma_ch4_bbody_counts']

    else:
        Ce = mmd['avhrr-ma_ch5_earth_counts']
        Cs = mmd['avhrr-ma_ch5_space_counts']
        Cict = mmd['avhrr-ma_ch5_bbody_counts']

    Tict = mmd['avhrr-ma_ict_temp'] # equivalent to mmd['avhrr-ma_orbital_temp']
    PRT1 = mmd['avhrr-ma_prt_1'][:,3,3]
    PRT2 = mmd['avhrr-ma_prt_2'][:,3,3]
    PRT3 = mmd['avhrr-ma_prt_3'][:,3,3]
    PRT4 = mmd['avhrr-ma_prt_4'][:,3,3]
    PRT = np.mean(np.vstack([PRT1, PRT2, PRT3, PRT4]).T, axis=1)

    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(PRT1, '.', markersize=0.5, label='prt_1')
    ax[0,0].plot(PRT2, '.', markersize=0.5, label='prt_2')
    ax[0,0].plot(PRT3, '.', markersize=0.5, label='prt_3')
    ax[0,0].plot(PRT4, '.', markersize=0.5, label='prt_4')
    ax[0,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,0].plot(Tict[:,3,3] - PRT, 'k.', markersize=0.5, label='Tict - <PRT>')
    ax[1,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[0,1].hist(PRT, bins=100, label='<PRT>')
    ax[0,1].hist(Tict[:,3,3],bins=100, label='Tict')
    ax[0,1].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,1].plot(Tict[:,3,3], PRT, 'k.',markersize=0.5)
    ax[1,1].set_aspect('equal', adjustable='box')
    ax[1,1].set_xlabel('Tict')
    ax[1,1].set_ylabel('<PRT>')    
    plt.tight_layout()
    plotfile = str(ch) + '_' + 'PRT.png'
    plt.savefig(plotfile)
    plt.close('all')

    #
    #   convert Tict --> Lict
    #

    Lict_LUT = con.bt2rad(Tict,channel,lut)
    Lict_CCI = con.bt2rad_cci(Tict,channel) 

    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(Lict_LUT[:,3,3], '.', markersize=0.5, label='Lict_LUT')
    ax[0,0].plot(Lict_CCI[:,3,3], '.', markersize=0.5, label='Lict_CCI')
    ax[0,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,0].plot(Lict_CCI[:,3,3] - Lict_LUT[:,3,3], 'k.', markersize=0.5, label='Lict_CCI - Lict_LUT')
    ax[1,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[0,1].hist(Lict_LUT[:,3,3], bins=100, label='Lict_LUT')
    ax[0,1].hist(Lict_CCI[:,3,3],bins=100, label='Lict_CCI')
    ax[0,1].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,1].plot(Lict_CCI[:,3,3], Lict_LUT[:,3,3], 'k.',markersize=0.5)
    ax[1,1].set_aspect('equal', adjustable='box')
    ax[1,1].set_xlabel('Lict_CCI')
    ax[1,1].set_ylabel('Lict_LUT')    
    plotfile = str(ch) + '_' + 'Lict.png'
    plt.tight_layout()
    plt.savefig(plotfile)
    plt.close('all')

    #
    #   convert Lict --> BTict
    #

    BTict_LUT = con.rad2bt(Lict_LUT,channel,lut)
    BTict_CCI = con.rad2bt_cci(Lict_LUT,channel)

    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(BTict_LUT[:,3,3], '.', markersize=0.5, label='BTict_LUT')
    ax[0,0].plot(BTict_CCI[:,3,3], '.', markersize=0.5, label='BTict_CCI')
    ax[0,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,0].plot(BTict_CCI[:,3,3] - BTict_LUT[:,3,3], 'k.', markersize=0.5, label='BTict_CCI - BTict_LUT')
    ax[1,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[0,1].hist(BTict_LUT[:,3,3], bins=100, label='BTict_LUT')
    ax[0,1].hist(BTict_CCI[:,3,3],bins=100, label='BTict_CCI')
    ax[0,1].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)    
    ax[1,1].plot(BTict_CCI[:,3,3], BTict_LUT[:,3,3], 'k.',markersize=0.5)
    ax[1,1].set_aspect('equal', adjustable='box')
    ax[1,1].set_xlabel('BTict_CCI')
    ax[1,1].set_ylabel('BTict_LUT')    
    plotfile = str(ch) + '_' + 'BTict.png'
    plt.tight_layout()
    plt.savefig(plotfile)
    plt.close('all')

    #
    # Test CCI conversion over full range of LUT
    #

    LUT_L = lut['L'][:,channel]
    LUT_BT = lut['BT'][:,channel]
    CCI_BT = con.rad2bt_cci(LUT_L,channel)

    fig, ax = plt.subplots()
    plt.plot(CCI_BT, CCI_BT-LUT_BT)
    plt.xlabel('BT (CCI) / $K$')
    plt.ylabel('BT difference (CCI-LUT) / $K$')
    plotfile = str(ch) + '_' + 'BT_CCI_v_LUT.png'
    plt.savefig(plotfile)
    plt.close('all')

    #
    #   convert Ce,Cs,Cict,Lict --> L
    #

    L_CCI = con.counts2rad_cci(channel,Ce,Cs,Cict,Lict_CCI)    

    #-------------------------------------------------------------------------------
    # Calculate radiance from counts and temperatures with measurement equation: HAR
    #-------------------------------------------------------------------------------

    #
    # Load: harmonisation coefficients
    #

    nsensor = 9
    if channel == 3:
        npar = 3
    else:
        npar = 4
    parameters = np.empty(shape=(npar*nsensor))

    har = xarray.open_dataset(har_file, decode_cf=True)

    parameter = har.parameter
    parameter_count = har.sensor_equation_parameter_count[1:].values
    parameter_pos = np.cumsum(har.sensor_equation_parameter_count).values
    for idx in range(len(parameter_count)):
        parameters[(idx*npar):(idx*npar+parameter_count[idx])] = parameter[(parameter_pos[idx]):(parameter_pos[idx+1])]
                     
    if channel == 3:
        a0 = parameters[(idx_ *npar)]
        a1 = parameters[(idx_ *npar)+1]
        a2 = parameters[(idx_ *npar)+2]
        a3 = 0.0
        a4 = 0.0
        if noT:
            a2 = 0.0
    else:
        a0 = parameters[(idx_ *npar)]
        a1 = parameters[(idx_ *npar)+1]
        a2 = parameters[(idx_ *npar)+2]
        a3 = parameters[(idx_ *npar)+3]
        a4 = 0.0
        if noT:
            a3 = 0.0

    #
    # Measurement equation correction term: Tinst
    #
    
    T_mean = np.mean(Tict[:,3,3])
    T_sdev = np.std(Tict[:,3,3])
    Tinst = (mmd['avhrr-ma_orbital_temperature'][:,3,3] - T_mean) / T_sdev

    #
    # Measurement equation correction term: WV (dummy for now)
    #

    WV = []

    Lict = Lict_CCI
    L_HAR = con.count2rad(Ce,Cs,Cict,Lict,Tinst,WV,channel,a0,a1,a2,a3,a4,noT)

    #
    # Convert radiance to BT
    #

    BT_CCI = con.rad2bt_cci(L_CCI,channel)[:,3,3]
    BT_HAR = con.rad2bt_cci(L_HAR,channel)[:,3,3]
    BT_MMD = BT_MMD[:,3,3]

    gd = BT_HAR > 0
    
    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(BT_CCI[gd],'.',markersize=0.5, label='BT(CCI)')
    ax[0,0].plot(BT_HAR[gd],'.',markersize=0.5, label='BT(HAR)')
    ax[0,0].plot(BT_MMD[gd],'.',markersize=0.5, label='BT(MMD)')
    ax[0,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)
    ax[0,0].set_ylim(200,320)
    ax[0,0].set_ylabel('BT / $K$')
    ax[1,0].plot(BT_MMD[gd], BT_MMD[gd]-BT_CCI[gd],'.',markersize=0.5, label='BT(MMD)-BT(CCI)')
    ax[1,0].plot(BT_MMD[gd], BT_MMD[gd]-BT_HAR[gd],'.',markersize=0.5, label='BT(MMD)-BT(HAR)')
    ax[1,0].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)
    ax[1,0].set_xlim(200,320)
    ax[1,0].set_ylim(-2,2)
    ax[1,0].set_xlabel('BT(MMD) / $K$')
    ax[1,0].set_ylabel('BT(MMD)-BT(x) / $K$')
    ax[0,1].hist((BT_MMD[gd]-BT_CCI[gd]),bins=100, label='BT(MMD)-BT(CCI)')  
    ax[0,1].hist((BT_MMD[gd]-BT_HAR[gd]),bins=100, label='BT(MMD)-BT(HAR)') 
    ax[0,1].set_xlim(-2,2)
    ax[0,1].set_xlabel('BT(MMD)-BT(x) / $K$')
    ax[0,1].set_ylabel('count (nbins=100)')
    ax[0,1].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)
    ax[1,1].plot(BT_MMD[gd], BT_CCI[gd], '.',markersize=0.5, label='BT(CCI)')
    ax[1,1].plot(BT_MMD[gd], BT_HAR[gd], '.',markersize=0.5, label='BT(HAR)')
    ax[1,1].legend(loc=1, markerscale=20, scatterpoints=5, fontsize=8)
    ax[1,1].set_xlim(200,320)
    ax[1,1].set_ylim(200,320)
    ax[1,1].set_xlabel('BT(MMD) / $K$')
    ax[1,1].set_ylabel('BT(x) / $K$')
    plt.tight_layout()
    plotfile = str(ch) + '_' + 'BT_v_BT_MMD.png'
    plt.savefig(plotfile)

    print('** END')




