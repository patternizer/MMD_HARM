import numpy as np
from netCDF4 import Dataset




def SST_time_anomaly_single(cloud,u,lat,month,hour):

    # SST_TIME_ANOMALY_SINGLE function
    # This function is for a single value anomaly calculation
    #
    # For array calculations, see SST_TIME_ANOMALY_ARRAY.
    #
    # Inputs required:
    # cloud - single value of cloud cover at time of comparison (0-1)
    # u - single value of wind speed (magnitude only) in m/s
    # lat - single value latitude of comparison
    # month - single number (1-12) denoting month of the year of comparison
    # hour - nearest hour (0-24) of the comparison



    if int(month) == 12 or int(month) > 2.5:
        seas = 0
    elif int(month) > 2.5 and int(month) < 5.5:
        seas = 1
    elif int(month) > 5.5 and int(month) < 8.5:
        seas = 2
    elif int(month) > 8.5 and int(month) < 11.5:
        seas = 3
    else:
       print "Error in month - not a number between 0-12"
       error = np.nan
       return error

    latbnd = np.arange(-90.0,91.0,10.0)
    
    for j in xrange(0,len(latbnd)):
        if latbnd[j] < lat and latbnd[j+1] > lat:
            lat_index = j
            break
    #### READ IN A0-A13 ##########
    data = Dataset('diurnal_sst_anomalies_from_drifting_buoys_v1.1.nc')
    coeffs = data.variables['fit_coeffs'][:]
    tod = data.variables['tod'][:]
    uavgdata = data.variables['u'][:]
    uavg = np.nanmean(uavgdata[lat_index,:,seas,:],axis=0)
    a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 = coeffs[lat_index,seas,:] 
    w = (2*np.pi)/24.0
    SSTanomalyday = np.zeros([24])
    i=0
    for t in tod:
        #print t, uavg[i], a2,w,u
        #print a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13
        #print np.shape(a0 + a1*u + a2*cloud + a3*t+ a4*u*t + a5*cloud*t + ((a6*(np.sin(w*t))) + a7*cloud*(np.sin(w*t)) + a8*(np.cos(w*t)) + a9*cloud*(np.cos(w*t)) + a10*(np.sin(2*w*t)) + a11*cloud*(np.sin(2*w*t)) + a12*(np.cos(2*w*t)) + a13*cloud*(np.cos(2*w*t)))*np.exp((-1.0*u)/uavg[i]))
        SSTanomalyday[i] = a0 + a1*u + a2*cloud + a3*t+ a4*u*t + a5*cloud*t + ((a6*(np.sin(w*t))) + a7*cloud*(np.sin(w*t)) + a8*(np.cos(w*t)) + a9*cloud*(np.cos(w*t)) + a10*(np.sin(2*w*t)) + a11*cloud*(np.sin(2*w*t)) + a12*(np.cos(2*w*t)) + a13*cloud*(np.cos(2*w*t)))*np.exp((-1.0*u)/uavg[i])

        i = i+1
    
    SST_anom = np.interp(hour,tod,SSTanomalyday)
    

    return SST_anom












def SST_time_anomaly_array(cloud,u,lat,month,hour):

    # SST_TIME_ANOMALY_ARRAY function
    # This function is for ARRAY value anomaly calculation
    # ALL ARRAYS MUST BE THE SAME SIZE WHEN ENTERED INTO THIS FUNCTION
    #
    # For SINGLE calculations, see SST_TIME_ANOMALY_SINGLE.
    #
    # Inputs required:
    # cloud - 1d array of cloud cover at time of comparison (0-1)
    # u - 1d array of wind speed (magnitude only) in m/s
    # lat - 1d array of latitude of comparison
    # month - 1d array of numbers(1-12) denoting month of the year of comparison
    # hour - 1d array of nearest hour (0-24) of the comparison


    comp_len = len(cloud)
    SST_anom = np.empty(comp_len)
    SST_anom[:] = np.nan
    data = Dataset('diurnal_sst_anomalies_from_drifting_buoys_v1.1.nc')
    coeffs = data.variables['fit_coeffs'][:]
    tod = data.variables['tod'][:]   
    uavgdata = data.variables['u'][:]  
    w = (2*np.pi)/24.0
    for comp in xrange(0,comp_len):
        if int(month[comp]) == 12 or int(month[comp]) > 2.5:
            seas = 0
        elif int(month[comp]) > 2.5 and int(month[comp]) < 5.5:
            seas = 1
        elif int(month[comp]) > 5.5 and int(month[comp]) < 8.5:
            seas = 2
        elif int(month[comp]) > 8.5 and int(month[comp]) < 11.5:
            seas = 3
        else:
           print "Error in month - not a number between 0-12"
           error = np.nan
           return error

        latbnd = np.arange(-90.0,91.0,10.0)

        for j in xrange(0,len(latbnd)):
            if latbnd[j] < lat[comp] and latbnd[j+1] > lat[comp]:
                lat_index = j
                break
        #### READ IN A0-A13 ##########        
        uavg = np.nanmean(uavgdata[lat_index,:,seas,:],axis=0)
        a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 = coeffs[lat_index,seas,:] 

        SSTanomalyday = np.zeros([24])
        i=0
        comp = int(comp)
        for t in tod:
            SSTanomalyday[i] = a0 + a1*u[comp] + a2*cloud[comp] + a3*t+ a4*u[comp]*t + a5*cloud[comp]*t + (((a6*(np.sin(w*t))) + a7*cloud[comp]*(np.sin(w*t)) + a8*(np.cos(w*t)) + a9*cloud[comp]*(np.cos(w*t)) + a10*(np.sin(2*w*t)) + a11*cloud[comp]*(np.sin(2*w*t)) + a12*(np.cos(2*w*t)) + a13*cloud[comp]*(np.cos(2*w*t)))*np.exp((-1*u[comp])/uavg[i]))

            i = i+1
    
        SST_anom[comp] = np.interp(hour[comp],tod,SSTanomalyday)
    

    return SST_anom



import datetime as dt


def secs_to_hour(seconds):

    # sec_to_hour function
    #
    # Turns either single value or 1d array into month,hours,seconds
    #
    # REQUIRES: seconds - either 1d array or single value of seconds since epoch
    #
    # RETURNS: month - single value int/1d array 1 - 12 of month
    #          hour - single value int/1d array of hour (0 - 23)
    #          minute - single value int/1d array of minute (0-59)

    if hasattr(seconds, "__len__"):
        #assume array
        array_len = len(seconds)
        month = np.zeros(array_len)
        hour = np.zeros(array_len)
        minute = np.zeros(array_len)
        i = 0
        for sec in seconds:
            month[i] = dt.datetime.fromtimestamp(sec).strftime('%m')
            hour[i] = dt.datetime.fromtimestamp(sec).strftime('%H')
            minute[i] = dt.datetime.fromtimestamp(sec).strftime('%M')
            i = i+1
    else:
       month = int(dt.datetime.fromtimestamp(seconds).strftime('%m'))
       hour = int(dt.datetime.fromtimestamp(seconds).strftime('%H'))
       minute = int(dt.datetime.fromtimestamp(seconds).strftime('%M'))

    return month,hour,minute





def SST_time_anomaly_comparison_array(cloud,u,lat,lon,month,hour,minute):

    # SST_TIME_ANOMALY_COMPARISON_ARRAY function
    # This function is for ARRAY value anomaly calculation
    # ALL ARRAYS MUST BE THE SAME SIZE WHEN ENTERED INTO THIS FUNCTION
    #
    #
    # Inputs required:
    # cloud - nx2 array of cloud cover at time of comparison (0-1)
    # u - nx2 array of wind speed (magnitude only) in m/s
    # lat - nx2 array of latitude of comparison
    # lon - nx2 array of longitude of comparison
    # month - nx2 array of numbers(1-12)denoting month of the year of comparison
    # hour - n x 2 array of hour (0-24) of the comparison
    # minute - n x 2 array of minute
    # 
    # In all cases, the 1st column should relate to one dataset for comparison
    # and the 2nd column, the other dataset for comparison
    #
    # 
    #
    # RETURNS: SST_anom - nx2 array of SST anomaly due to time
    #          SST_diff - n sized array of the difference between SST_anom for
    #                     the same comparison site.

    comp_len = len(cloud)
    data = Dataset('diurnal_sst_anomalies_from_drifting_buoys_v1.1.nc')
    coeffs = data.variables['fit_coeffs'][:]
    tod = data.variables['tod'][:]   
    uavgdata = data.variables['u'][:]
    SST_anom = np.empty([comp_len,2])  
    SST_anom[:,:] = np.nan

    w = (2*np.pi)/24.0

    for col in xrange(0,2):
        for comp in xrange(0,comp_len):
            if int(month[comp,col]) == 12 or int(month[comp,col]) < 2.5:
                seas = 0
            elif int(month[comp,col]) > 2.5 and int(month[comp,col]) < 5.5:
                seas = 1
            elif int(month[comp,col]) > 5.5 and int(month[comp,col]) < 8.5:
                seas = 2
            elif int(month[comp,col]) > 8.5 and int(month[comp,col]) < 11.5:
                seas = 3
            else:
               print "Error in month - not a number between 0-12"
               error = np.nan
               return error

            latbnd = np.arange(-90.0,91.0,10.0)

            for j in xrange(0,len(latbnd)):
                if latbnd[j] < lat[comp,col] and latbnd[j+1] > lat[comp,col]:
                    lat_index = j
                    break
            #### READ IN A0-A13 ##########        
            uavg = np.nanmean(uavgdata[lat_index,:,seas,:],axis=0)
            a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13 = coeffs[lat_index,seas,:] 

            SSTanomalyday = np.zeros([24])
            i=0
            comp = int(comp)
            for t in tod:
                SSTanomalyday[i] = a0 + a1*u[comp,col] + a2*cloud[comp,col] + a3*t+ a4*u[comp,col]*t + a5*cloud[comp,col]*t + (((a6*(np.sin(w*t))) + a7*cloud[comp,col]*(np.sin(w*t)) + a8*(np.cos(w*t)) + a9*cloud[comp,col]*(np.cos(w*t)) + a10*(np.sin(2*w*t)) + a11*cloud[comp,col]*(np.sin(2*w*t)) + a12*(np.cos(2*w*t)) + a13*cloud[comp,col]*(np.cos(2*w*t)))*np.exp((-1*u[comp,col])/uavg[i]))

                i = i+1
            min_decimal = float(minute[comp,col])/60.0
            UTC_time = hour[comp,col]+min_decimal
            local_time = UTC_time + (lon[comp,col]/15.0)
            SST_anom[comp,col] = np.interp((local_time),tod,SSTanomalyday,period=24.0)
    np.shape(SST_anom)
    SST_diff = (SST_anom[:,0] - SST_anom[:,1])
    return SST_anom, SST_diff




