#!/usr/bin/env python

# Code to quickly plot L1C orbits from FIDUCEO Easy FCDR:
# ======================================================
# Version 0.4
# 15 June, 2019
# michael.taylor AT reading DOT ac DOT uk
# patternizer.github.io
# ======================================================

from  optparse import OptionParser
import numpy as np
import numpy.ma as ma
import xarray
import matplotlib.pyplot as plt; plt.close("all")
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#####################################################################################
# plot_orbit_var(lat, lon, var, vmin, vmax, projection, filestr, titlestr, varstr)
#####################################################################################                              
#
# REQUIRES:                               
# 
# lat: np.array
# lon: np.array
# var: np.array
# vmin: e.g. 270
# vmax: e.g. 305
#
# projection: string - one of:
# ['platecarree','molleweide','robinson','equalearth','geostationary',
# 'goodehomolosine','europp','northpolarstereo','southpolarstereo','lambertconformal']
#
# filestr: string (including file extension) e.g. 'plot_l1c.png'
# titlestr: string e.g. 'L1C Orbit'
# varstr: e.g. 'ch4'
#
#####################################################################################

def plot_orbit_var(lat, lon, var, vmin, vmax, projection, filestr, titlestr, varstr):

    x = lon[::10,::10]
    y = lat[::10,::10]
    z = var[::10,::10]
        
    fig  = plt.figure()
    if projection == 'platecarree':
        p = ccrs.PlateCarree(central_longitude=0)
        threshold = 0
    if projection == 'mollweide':
        p = ccrs.Mollweide(central_longitude=0)
        threshold = 1e6
    if projection == 'robinson':
        p = ccrs.Robinson(central_longitude=0)
        threshold = 0
    if projection == 'equalearth':
        p = ccrs.EqualEarth(central_longitude=0)
        threshold = 0
    if projection == 'geostationary':
        p = ccrs.Geostationary(central_longitude=0)
        threshold = 0
    if projection == 'goodehomolosine':
        p = ccrs.InterruptedGoodeHomolosine(central_longitude=0)
        threshold = 0
    if projection == 'europp':
        p = ccrs.EuroPP()
        threshold = 0
    if projection == 'northpolarstereo':
        p = ccrs.NorthPolarStereo()
        threshold = 0
    if projection == 'southpolarstereo':
        p = ccrs.SouthPolarStereo()
        threshold = 0
    if projection == 'lambertconformal':
        p = ccrs.LambertConformal(central_longitude=0)
        threshold = 0

    ax = plt.axes(projection=p)    
    ax.stock_img()
    ax.coastlines(resolution='50m')
    ax.gridlines()
    colormap = 'gnuplot2'

    g = ccrs.Geodetic()
    trans = ax.projection.transform_points(g, x, y)
    x0 = trans[:,:,0]
    x1 = trans[:,:,1]
    if projection == 'platecarree':
        ax.set_extent([-180, 180, -90, 90], crs=p)
        gl = ax.gridlines(crs=p, draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='-')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.ylines = True
        gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
        gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        for mask in (x0>=threshold,x0<=threshold):
            im = ax.pcolor(ma.masked_where(mask, x), ma.masked_where(mask, y), ma.masked_where(mask, z), vmin=vmin, vmax=vmax, transform=ax.projection, cmap=colormap)
    else:
        for mask in (x0>=threshold,x0<=threshold):
            im = ax.pcolor(ma.masked_where(mask, x0), ma.masked_where(mask, x1), ma.masked_where(mask, z), vmin=vmin, vmax=vmax, transform=ax.projection, cmap=colormap)
    im.set_clim(vmin,vmax)
    cb = plt.colorbar(im, orientation="horizontal", extend='both', label=varstr)
    plt.title(titlestr)
    plt.savefig(filestr)
    plt.close('all')

    return


