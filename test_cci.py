#!/usr/bin/env python

# Compare L1C MMD and ME brightness temperatures
# ======================================================
# Version 0.4
# 18 June, 2019
# michael.taylor AT reading DOT ac DOT uk
# patternizer.github.io
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

    ################################################################
    # RUN PARAMETERS: 
    ################################################################

    ch = 12
    idx = 7      # MTA (see avhrr_sat)
    cci_file = 'mta_cci.nc'
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
    # senesor list and indexing
    #

    avhrr_sat = [b'N12',b'N14',b'N15',b'N16',b'N17',b'N18',b'N19',b'MTA',b'MTB']

    # NB: Harmonisation coefficients provided by Ralf Quast are for MTA --> N11:
    # RQ:          MTA,   N19,   N18,   N17,   N16,   N15,  N14,    N12,   N11
    # index         0       1      2      3      4      5      6     7       8
    # --> new index map

    idx_ = 7 - idx

    #
    # Load: lut
    #

    lut = con.read_in_LUT(avhrr_sat[idx])

    #
    # Load: ESA_CCI L1C (v2.10.1) temperatures
    #

    cci = xarray.open_dataset(cci_file, decode_times=False)

    if channel == 3: BT_CCI = cci['ch3b'][0,:,:]
    elif channel == 4: BT_CCI = cci['ch4'][0,:,:]
    else: BT_CCI = cci['ch5'][0,:,:]

    ##############################
    # v1.5 CCI-like coefficients #
    ##############################
     
    # [3.7]:   L = a1 + ( (Lict * (0.985140 + a2) ) / (Cict - Cs) ) * (Ce - Cs) + a3 * Tinst
    # [11,12]: L = a1 + ( (Lict * (0.985140 + a2) ) / (Cict - Cs) + a3 * (Ce - Cict) ) * (Ce - Cs) + a4 * Tinst
    # -->      L = a1 + ( (Lict * (0.985140 + a2) ) / (Cict - Cs) ) * (Ce - Cs) + a4 * Tinst 
    #                 + a3 * (Ce - Cict) * (Ce - Cs)

    # Metop-A non-linear coefficients etc.:
    # --------------------------------------------------------------------------------------
    # https://github.com/surftemp/gbcs/blob/master/src/AVHRR/AVHRR_Calibration.f90#L464:L470

    # CASE(-1)
    # year_launch= 2006.77
    # ch1_gain_l = [0.056, 0.906, -0.024]
    # ch1_gain_h = [0.167, 0.906, -0.024]
    # ch2_gain_l = [0.066, 0.814, 0.025]
    # ch2_gain_h = [0.199, 0.814, 0.025]
    # ch3_gain_l = [0.031, 2.020, -0.115]
    # ch3_gain_h = [0.220, 2.020, -0.115]
    # coef%vis_switch = [501.00, 500.00, 502.00]

    # coef%Central_Wavenumber = [2687.0392, 927.27630, 837.80762]
    # coef%Const_offset   = [-2.0653147, -0.56503332, -0.38472766]
    # coef%Const_slope    = [1.0034418, 1.0015090, 1.0011264]
    # coef%NSpace         = [0., -4.98, -3.40]
    # coef%nonlinear(:,1) = [0., 0., 0.]
    # coef%nonlinear(:,2) = [5.44,-0.10152,0.00046964]
    # coef%nonlinear(:,3) = [3.84,-0.06249,0.00025239]
 
    # Calculate the TIR gain (assumes you've filtered and smoothed space/bbody etc):
    # ---------------------------------------------------------------------------------------
    # https://github.com/surftemp/gbcs/blob/master/src/AVHRR/AVHRR_Calibration.f90#L990:L1023

    # PURE SUBROUTINE Calculate_TIR_Gain(caldata)
    # REAL(GbcsReal) :: cwn, pc1, pc2, tstar, r_ict, diff

    # DO ichan=1,3
    #   ! Temperature to radiance conversion coefficients
    #   cwn = caldata%coef%Central_Wavenumber(ichan)
    #   pc1 = Planck_C1 * cwn**3!*cwn*cwn
    #   pc2 = Planck_C2 * cwn
    #   DO iscan=1, caldata%nlines
    #     IF (caldata%t_ict(iscan) < 0                .OR. &
    #         caldata%smooth_space(iscan,ichan+3) < 0 .OR. &
    #         caldata%smooth_bbody(iscan,ichan+3) < 0) CYCLE
    #     tstar = (caldata%t_ict(iscan) - caldata%coef%Const_Offset(ichan)) / caldata%coef%Const_Slope(ichan)
    #     r_ict = pc1 / (EXP(pc2/tstar) - 1)
    #     ! bbody and space are both allocated for six channels
    #     diff = caldata%smooth_space(iscan,ichan+3) - caldata%smooth_bbody(iscan,ichan+3)
    #     IF (diff > 10) THEN
    #       caldata%gain(iscan,ichan) = REAL((r_ict - caldata%coef%NSPace(ichan)) / diff)
    #     END IF
    #   END DO
    # END DO

    # Calculate the TIR calibration coefficients:
    # -----------------------------------------------------------------------------------------
    # https://github.com/surftemp/gbcs/blob/master/src/AVHRR/AVHRR_Calibration.f90#L1061::L1092

    # SUBROUTINE Calculate_TIR_Cal(caldata)
    # REAL(GbcsDble) :: Nspace
    # REAL(GbcsDble), DIMENSION(3) :: nonlin
    # REAL(GbcsDble) :: a0, a1, a2, s, g

    # DO ichan=1,3
    #   NSpace = caldata%coef%NSpace(ichan)
    #   nonlin = caldata%coef%nonlinear(:,ichan)
    #   a0 = nonlin(1) + (1 + nonlin(2) + nonlin(3)*Nspace)*Nspace
    #   a1 = (1 + nonlin(2) + 2*nonlin(3)*Nspace)
    #   a2 = nonlin(3)
    #   DO iscan=1, caldata%nlines
    #     s = caldata%smooth_space(iscan,ichan+3)
    #     g = caldata%gain(iscan,ichan)
    #     IF (s < 0 .OR. g < 0) CYCLE
    #     caldata%tir_cal(1,iscan,ichan) = REAL(a0 + a1*s*g + a2*s*s*g*g)
    #     caldata%tir_cal(2,iscan,ichan) = REAL(-a1*g - 2*a2*s*g*g)
    #     caldata%tir_cal(3,iscan,ichan) = REAL(a2*g*g)
    #   END DO
    # END DO

    # Convert earth counts to radiance (and then to BT):
    # ----------------------------------------------------------------------------------------
    # https://github.com/surftemp/gbcs/blob/master/src/AVHRR/AVHRR_Calibration.f90#L1541:L1582

    # CASE(Channel_3B, Channel_4, Channel_5)
    #  ! -- IR channel
    #  tir_cal = cal%tir_cal(:,scan,chan-3)
    #  IF (tir_cal(1) < 0) THEN
    #    ! Bad calibration
    #    out = NAN_R
    #    err = NAN_R
    #    RETURN
    #  END IF

    # cwn   = cal%coef%Central_Wavenumber(chan-3)
    # offset= cal%coef%Const_offset(chan-3)
    # slope = cal%coef%Const_slope(chan-3)
    # coef1 = Planck_C1 * cwn*cwn*cwn
    # coef2 = Planck_C2 * cwn

    # DO elem=1, cal%nelems
    #   counts = avhrrdata(scan)%earth_counts(chan-1,elem)
    #   rads = tir_cal(1) + tir_cal(2)*counts + tir_cal(3)*counts*counts

    #   IF (cal%apply_radbias) THEN
    #     ! Used for Walton Tinstr and scene temperature bias adjustments
    #     rads = cal%radiance_bias(1,chan-3) + cal%radiance_bias(2,chan-3) * rads
    #   END IF

    #   IF (rads <= 0) THEN
    #     out(elem) = NAN_R
    #     err(elem) = NAN_R
    #     CYCLE
    #   END IF

    #   ! -- Convert radiometric noise from counts to radiance
    #   noise = (tir_cal(2) + 2*tir_cal(3)*counts) * noise_counts
    #   ! -- Add in systematic component (as we only have one error component)
    #   noise = SQRT(noise**2 + cal%uncert_sys(chan)**2)

    #   out(elem) = Calc_BT(rads, coef1, coef2, offset, slope)
    #   ! -- Convert delta-rads to delta-temperature
    #   tstar = (out(elem) - offset) / slope
    #   noise = noise * coef1 * tstar * tstar / (coef2 * rads * (rads + coef1))
    #   err(elem) = noise
    # END DO      

    #
    # Load: harmonisation coefficients
    #

    ds = xarray.open_dataset(harm_file, decode_cf=True)

    # OLD & NEW FORMATS:
    # -----------------

    # OLD: <xarray.DataArray 'parameter_sensors' (n: 36)>
    # array([b'm02',b'm02',b'm02',b'm02',
    # b'n19',b'n19',b'n19',b'n19',
    # b'n18',b'n18',b'n18',b'n18',
    # b'n17',b'n17',b'n17',b'n17',
    # b'n16',b'n16',b'n16',b'n16',
    # b'n15',b'n15',b'n15',b'n15',
    # b'n14',b'n14',b'n14',b'n14',
    # b'n12',b'n12',b'n12',b'n12',
    # b'n11',b'n11',b'n11',b'n11'])

    # NEW: <xarray.DataArray 'parameter_name' (n: 34)>
    # array([b'm02',b'm02',b'm02',
    # b'n19',b'n19',b'n19',
    # b'n18',b'n18',b'n18',b'n18',
    # b'n17',b'n17',b'n17',b'n17',
    # b'n16',b'n16',b'n16',b'n16',
    # b'n15',b'n15',b'n15',b'n15',
    # b'n14',b'n14',b'n14',b'n14',
    # b'n12',b'n12',b'n12',b'n12',
    # b'n11',b'n11',b'n11',b'n11'])

    # OLD: N/A

    # NEW: <xarray.DataArray 'sensor_name' (n_sensor: 10)>
    # array([b'aatsr',b'm02',b'n19',b'n18',b'n17',b'n16',b'n15',b'n14',b'n12',b'n11'])

    # OLD: <xarray.DataArray 'parameter_count' (s: 9)>
    # array([4, 4, 4, 4, 4, 4, 4, 4, 4], dtype=int32)

    # NEW: <xarray.DataArray 'sensor_equation_parameter_count' (n_sensor: 10)>
    # array([0, 3, 3, 4, 4, 4, 4, 4, 4, 4], dtype=int32)    

    # OLD: <xarray.DataArray 'parameter' (n: 36)>
    # array([ 9.530414e-01,  8.612669e-03,  2.694954e-06, -2.812751e-03, 
    # 3.153889e-01,  3.026021e-03,  5.301578e-06, -4.093177e-02,
    # 7.927444e-01,  3.160753e-03,  2.611367e-06,  3.073914e-01,
    # 9.591988e-01,  1.007047e-02,  3.220955e-06, -4.815025e-02,
    # 1.062143e+00,  2.624904e-03, -6.036055e-07,  1.793575e-01,
    # 1.047737e+00,  2.167977e-03, -8.602448e-06,  5.509021e-02,
    # 1.107450e+00,  2.227697e-02,  9.925325e-06,  1.676749e-01,
    # 8.568336e-01,  3.673877e-02, -3.724639e-06,  2.174961e-01,
    # 1.919452e+00,  5.799905e-02,  3.035878e-05,  2.117777e-01])

    # NEW: <xarray.DataArray 'parameter' (n: 34)>
    # array([ 1.042599e+00,  9.176152e-03,  4.827715e-06,  5.051670e-01,
    #         3.914914e-03,  9.023554e-06,  9.372754e-01,  4.200608e-03,
    #         5.571235e-06,  3.344755e-01,  1.064974e+00,  1.096942e-02,
    #         5.546465e-06, -1.587046e-02,  1.149384e+00,  4.419415e-03,
    #         1.911674e-06,  1.974459e-01,  1.244548e+00,  5.234328e-03,
    #        -3.590035e-06,  7.734781e-02,  1.668203e+00,  2.612718e-02,
    #         1.965677e-05,  2.243123e-01,  1.487629e+00,  4.219526e-02,
    #         7.727572e-06,  2.852593e-01,  2.562915e+00,  6.439369e-02,
    #         4.307746e-05,  3.219609e-01])

    # OLD: <xarray.DataArray 'parameter_add_offset' (n: 36)>    
    # array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    # 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

    # NEW: <xarray.DataArray 'parameter_add_offset' (n: 34)>
    # array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    # 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


    # OLD: <xarray.DataArray 'parameter_scale_factor' (n: 36)>
    # array([1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00])

    # NEW: <xarray.DataArray 'parameter_scale_factor' (n: 34)>
    # array([1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e-03, 1.e-06, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00, 1.e+00, 1.e-03, 1.e-06, 1.e+00, 
    # 1.e+00, 1.e-03, 1.e-06, 1.e+00])
       
    # OLD: <xarray.DataArray 'measurement_equation_id' (s: 9)>
    # array([102, 102, 102, 102, 102, 102, 102, 102, 102], dtype=int32)

    # NEW: <xarray.DataArray 'sensor_equation_id' (n_sensor: 10)>
    # array([  0, 101, 101, 102, 102, 102, 102, 102, 102, 102], dtype=int32)

    # OLD: <xarray.DataArray 'measurement_equation_constant_count' (s: 9)>
    # array([2, 2, 2, 2, 2, 2, 2, 2, 2], dtype=int32)

    # NEW: <xarray.DataArray 'sensor_equation_config_count' (n_sensor: 10)>
    # array([0, 0, 0, 2, 2, 2, 2, 2, 2, 2], dtype=int32)
    
    # OLD: <xarray.DataArray 'measurement_equation_constant' (c: 18)>
    # array([2.861258e+02, 4.908800e-02, 2.877546e+02, 1.176810e-01, 
    # 2.882198e+02, 6.076970e-01, 2.881066e+02, 1.607656e+00, 
    # 2.926722e+02, 3.805704e+00, 2.947586e+02, 2.804361e+00, 
    # 2.886376e+02, 1.053762e+00, 2.903271e+02, 2.120666e+00, 
    # 2.904022e+02, 3.694937e+00])

    # NEW: <xarray.DataArray 'sensor_equation_config' (n_config: 14)>
    # array([288.219774,   0.607697, 288.10663 ,   1.607656, 292.672201,   3.805704,
    #        294.758564,   2.804361, 288.637636,   1.053762, 290.327113,   2.120666,
    #        290.402168,   3.694937])

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

    elif channel == 4:
        Ce = mmd['avhrr-ma_ch4_earth_counts']
        Cs = mmd['avhrr-ma_ch4_space_counts']
        Cict = mmd['avhrr-ma_ch4_bbody_counts']

    else:
        Ce = mmd['avhrr-ma_ch5_earth_counts']
        Cs = mmd['avhrr-ma_ch5_space_counts']
        Cict = mmd['avhrr-ma_ch5_bbody_counts']

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

    #===============================================

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
    plotfile = str(ch) + '_' + 'BT_v_BT_MMD.png'
    plt.savefig(plotfile)

    print('** END')


