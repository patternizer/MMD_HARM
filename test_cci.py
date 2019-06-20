#!/usr/bin/env python

# Compare L1C MMD and ME brightness temperatures
# ======================================================
# Version 0.6
# 20 June, 2019
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
    l1b_file = 'mta_l1b.nc'
    l1c_file = 'mta_l1c.nc'
    cci_file = 'mta_cci.nc'
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

    avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']

    # NB: Harmonisation coefficients provided by Ralf Quast:

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
    # Load: L1b orbit counts and temperatures
    #

    l1b = xarray.open_dataset(l1b_file, decode_times=False)

    #
    # Load: L1c orbit counts and temperatures
    #

    l1c = xarray.open_dataset(l1c_file, decode_times=False)

    #
    # Load: ESA_CCI L1C (v2.10.1) counts and temperatures
    #

    cci = xarray.open_dataset(cci_file, decode_times=False)

    #
    # Load: mmd orbit counts and temperatures
    #

    mmd = xarray.open_dataset(mmd_file, decode_times=False)

    if channel == 3: 
        BT_L1B = l1b['ch3b']             # (12348, 409)
        BT_L1C = l1c['Ch3b']             # (14062, 409)
        BT_CCI = cci['ch3b'][0,:,:]      # (12233, 409)
        BT_MMD = mmd['avhrr-ma_ch3b']    # (55604, 7, 7)
    elif channel == 4: 
        BT_L1B = l1b['ch4']
        BT_L1C = l1c['Ch4']
        BT_CCI = cci['ch4'][0,:,:]
        BT_MMD = mmd['avhrr-ma_ch4']
    else: 
        BT_L1B = l1b['ch5']
        BT_L1C = l1c['Ch5']
        BT_CCI = cci['ch5'][0,:,:]
        BT_MMD = mmd['avhrr-ma_ch5']

    ##############################
    # v1.5 CCI-like coefficients #
    ##############################
     
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
#    PRT = ((PRT1 + PRT2 + PRT3 + PRT4) / 4.)
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

    # METHOD 1: LUT

    Lict_LUT = con.bt2rad(Tict,channel,lut)

    # METHOD 2: CCI

    offset = np.append( np.zeros(3), np.array([-2.0653147, -0.56503332, -0.38472766]) )      # Aval
    slope = np.append( np.zeros(3), np.array([1.0034418, 1.0015090, 1.0011264]) )            # Bval
#   Central_Wavenumber = np.append( np.zeros(3), np.array([1687.0392, 927.2763, 837.80762])) # NuC (in CCI)
    Central_Wavenumber = np.append( np.zeros(3), np.array([2687.0392, 927.2763, 837.80762])) # NuC (in Jon's L1B)
    tstar = (Tict - offset[channel])/slope[channel]
    Planck_C1 = 0.00001191042722
    Planck_C2 = 1.4387752
    coef1 = Planck_C1 * Central_Wavenumber[channel]**3
    coef2 = Planck_C2 * Central_Wavenumber[channel]
    Lict_CCI = (coef1 / ( np.exp(coef2/tstar) - 1.0)) 

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

    BTict_LUT = con.rad2bt(Lict_LUT,channel,lut)
    BTict_CCI = con.rad2bt(Lict_CCI,channel,lut)

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

    #-------------------------------------------------------------------------------
    # Calculate radiance from counts and temperatures with measurement equation: CCI
    #-------------------------------------------------------------------------------

    Lict = Lict_CCI

    Nspace = np.append( np.zeros(3), np.array([0, -4.98, -3.4]) )
    nonlinear = np.vstack([np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3), np.array([5.44, -0.10152, 0.00046964]), np.array([3.84, -0.06249, 0.00025239])]).T
    a0 = nonlinear[0,channel] + (1. + nonlinear[1,channel] + nonlinear[2,channel]*Nspace[channel]) * Nspace[channel]
    a1 = 1. + nonlinear[1,channel] + 2. * nonlinear[2,channel] * Nspace[channel]
    a2 = nonlinear[2,channel]
    s = Cs
    g = (Lict - Nspace[channel]) / (Cs - Cict)
    term1 = a0 + a1*s*g + a2*s*s*g*g
    term2 = -1.*a1*g - 2.*a2*s*g*g
    term3 = a2*g*g
    L_CCI = (term1 + term2*Ce + term3*Ce*Ce)

    #-------------------------------------------------------------------------------
    # Calculate radiance from counts and temperatures with measurement equation: HAR
    #-------------------------------------------------------------------------------

    #
    # Load: harmonisation coefficients
    #

    har = xarray.open_dataset(har_file, decode_cf=True)

    par = np.cumsum(har.sensor_equation_parameter_count)
#    A = har.parameter(par[idx_ -1]:par[idx_]) 
    if channel == 37:
        a0 = har.parameter[(idx_ *3)].values
        a1 = har.parameter[(idx_ *3)+1].values
        a2 = har.parameter[(idx_ *3)+2].values
        a3 = 0.0
        a4 = 0.0
        if noT:
            a2 = 0.0
    else:
        a0 = har.parameter[(idx_ *4)].values
        a1 = har.parameter[(idx_ *4)+1].values
        a2 = har.parameter[(idx_ *4)+2].values
        a3 = har.parameter[(idx_ *4)+3].values
        a4 = 0.0
        if noT:
            a3 = 0.0

    #
    # Measurement equation correction term: Instrumental temperature
    #
    
#    T_mean = np.array([ 290.327113, 288.637636, 294.758564, 292.672201, 288.106630, 288.219774, 287.754638, 286.125823, np.nan])
#    T_sdev = np.array([ 2.120666, 1.053762, 2.804361, 3.805704, 1.607656, 0.607697, 0.117681, 0.049088, np.nan])

    T_mean = np.mean(Tict[:,3,3])
    T_sdev = np.std(Tict[:,3,3])
    Tinst = (mmd['avhrr-ma_orbital_temperature'] - T_mean[idx]) / T_sdev[idx]

    #
    # Measurement equation correction term: Water vapour (dummy for now)
    #

    WV = 0.0 * Tinst

    L_HAR = con.count2rad(Ce,Cs,Cict,Lict,Tinst,WV,channel,a0,a1,a2,a3,a4,noT) / 100.

    #
    # Convert radiance to BT
    #

    BT_LUT = con.rad2bt(L_LUT,channel,lut)
    BT_CCI = con.rad2bt(L_CCI,channel,lut)
    BT_HAR = con.rad2bt(L_HAR,channel,lut)
    
    #===============================================

    L_LUT = L_LUT[:,3,3]
    L_CCI = L_CCI[:,3,3]
    L_HAR = L_HAR[:,3,3]
    BT_CCI = BT_CCI[:,3,3]
    BT_HAR = BT_HAR[:,3,3]
    BT_MMD = BT_MMD[:,3,3]

    gd = BT_HAR > 0
    
    fig, ax = plt.subplots(2, 2) 
    ax[0,0].plot(BT_CCI[gd],'.',markersize=0.5, label='BT (from CCI)')
    ax[0,0].plot(BT_HAR[gd],'.',markersize=0.5, label='BT (from HAR)')
    ax[0,0].plot(BT_MMD[gd],'.',markersize=0.5, label='BT (from MMD)')
    ax[0,0].legend(markerscale=20, scatterpoints=5, fontsize=8)
    ax[0,0].set_ylim(200,320)
    ax[0,0].set_xlabel('MM')
    ax[0,0].set_ylabel('BT')
    ax[1,0].plot(BT_MMD[gd], BT_MMD[gd]-BT_CCI[gd],'.',markersize=0.5)
    ax[1,0].plot(BT_MMD[gd], BT_MMD[gd]-BT_HAR[gd],'.',markersize=0.5)
    ax[1,0].set_xlim(200,320)
    ax[1,0].set_ylim(-1,1)
    ax[1,0].set_xlabel('BT_MMD')
    ax[1,0].set_ylabel('BT(MMD)-BT(x)')
    ax[0,1].hist((BT_MMD[gd]-BT_CCI[gd]),bins=100)
    ax[0,1].hist((BT_MMD[gd]-BT_HAR[gd]),bins=100)
    ax[0,1].set_xlabel('BT(MMD)-BT(x)')
    ax[1,1].set_ylabel('count (nbins=100)')
    ax[1,1].plot(BT_MMD[gd], BT_CCI[gd], '.',markersize=0.5)
    ax[1,1].plot(BT_MMD[gd], BT_HAR[gd], '.',markersize=0.5)
    ax[1,1].set_xlim(200,320)
    ax[1,1].set_ylim(200,320)
    ax[1,1].set_xlabel('BT_MMD')
    ax[1,1].set_ylabel('BT(x)')
    plt.tight_layout()
    plotfile = str(ch) + '_' + 'BT_v_BT_MMD.png'
    plt.savefig(plotfile)

    print('** END')




