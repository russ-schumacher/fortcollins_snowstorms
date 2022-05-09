#!/usr/bin/env python
# coding: utf-8

# Code to make a series of ERA5 maps for Fort Collins snowstorms (or really any specified weather events)
# ----------
# 
# ### reads in a text file with a list of storms, then goes on to remotely access then plot various fields
# 
# #### updated to work with MetPy V1.X, spring 2022
# 

# #### Imports

# In[1]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import ctables
import numpy as np
import scipy.ndimage as ndimage
import cdsapi
import fsspec
import xarray as xr
import matplotlib.colors as mcolors
import os
import pandas as pd

from palettable.colorbrewer.diverging import BrBG_11


# In[2]:


## get our list of cases to work from
#case_date = '20191125'
case_data_all = pd.read_csv('snowstorms_12inch_fortcollins_post1979.csv')
#case_data_all = pd.read_csv('snowstorms_12inch_fortcollins_update.csv')

num_cases = case_data_all['start_date'].shape[0]
print('working on '+str(num_cases)+' cases')

### specify frequency of plots (hourly, every 6 hours, etc.)
mapfreq = 6  ## in hours

### specify number of days before and after event to plot
numdays = 4  ### e.g., if this is 2, plot 2 days before to 2 days after event


# In[3]:


## which fields to plot?
plot_250wsp = True
plot_500vort = True
plot_700wsp = True
plot_850temp = True
plot_mslp_thk = True
plot_mslp_precip = True
plot_700qvect = True
plot_pwat = True

plot_dt_theta = True

plot_sfc_dwpt = False
plot_mucape = False
plot_sfc_rh = False


# In[4]:


## set some plt defaults
plt.rcParams['savefig.facecolor'] = "white"
plt.rcParams['savefig.transparent'] = False
plt.rcParams['savefig.bbox'] = "tight"
plt.rcParams['savefig.dpi'] = 100

# Set Projection of Data
datacrs = ccrs.PlateCarree()

# Set Projection of Plots
plotcrs = ccrs.LambertConformal(central_latitude=45, central_longitude=-100)

## endash:
enDash = u'\u2013'


# In[5]:


## loop over cases
for cc in range(0,num_cases):
#for cc in range(19,num_cases):
#for cc in range(1,3):
#cc=0
    case_date = case_data_all['start_date'][cc]

    print(' ')
    print('working on '+str(case_date))

    #Data Acquisition

    case_date = pd.to_datetime(case_date, format="%Y-%m-%d")
    os.system("mkdir -p "+case_date.strftime("%Y%m%d"))

    plots_start_day = case_date - pd.Timedelta(days=numdays)
    plots_end_day = case_date + pd.Timedelta(days=numdays) + pd.Timedelta(hours=23)

    ### range of dates we want to plot
    dates = pd.date_range(start=plots_start_day, end=plots_end_day, freq=str(mapfreq)+"H")

    ### and for precip, need to start even a day earlier than this
    dates_precip = pd.date_range(start=(plots_start_day - pd.Timedelta(days=1)),end=plots_end_day)
    
    ### get the data from ECMWF API
    #### instead of downloading it, read it directly, based on: https://nbviewer.jupyter.org/github/pangeo-data/pangeo-era5/blob/master/lazy_loading_expts.ipynb

    print("getting pressure level data")

    c1 = cdsapi.Client()
    name = 'reanalysis-era5-pressure-levels'
    request = {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [
                'geopotential', 'relative_humidity', 'temperature',
                'u_component_of_wind', 'v_component_of_wind',
            ],
            'pressure_level': [
                '250', '500', '700',
                '850', '925', '1000',
            ],
            ### use the 'dates' array above to build the list of dates we want
            'date': dates[0].strftime("%Y-%m-%d")+"/"+dates[-1].strftime("%Y-%m-%d"),
            'time': dates.strftime("%H:%M").unique().tolist(),
            'area'          : [75.,185, 15., 320.] #N,W,S,E
    }
    r1 = c1.retrieve(name, request, None)

    with fsspec.open(r1.location) as f1:
        iso_vars = xr.open_dataset(f1, engine='scipy')
    #    print(iso_vars)
    
    iso_vars.to_netcdf("iso_vars_"+case_date.strftime("%Y%m%d")+".nc")

    print('getting single-level data')

    ## also get some single-level data
    ## get the data from ECMWF API
    c2 = cdsapi.Client()
    name = 'reanalysis-era5-single-levels'
    request = {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': ['mean_sea_level_pressure','total_column_water_vapour','2m_temperature',
                            '2m_dewpoint_temperature','10m_u_component_of_wind','10m_v_component_of_wind',
                            'convective_available_potential_energy'],
            'date': dates[0].strftime("%Y-%m-%d")+"/"+dates[-1].strftime("%Y-%m-%d"),
            'time': dates.strftime("%H:%M").unique().tolist(),
            'area'          : [75.,185, 15., 320.] #N,W,S,E
    }
    r2 = c2.retrieve(name, request, None)

    with fsspec.open(r2.location) as f2:
        sl_vars = xr.open_dataset(f2, engine='scipy')
        #print(sl_vars)
        
    sl_vars.to_netcdf("sl_vars_"+case_date.strftime("%Y%m%d")+".nc")

    ## now get precip accumulation from EC

    ## what's one day earlier because we need to download that too?
    if plot_mslp_precip:

        print('getting precip data')

        c3 = cdsapi.Client()
        name = 'reanalysis-era5-single-levels'
        request = {
                'variable'      : 'total_precipitation',
                'product_type': 'reanalysis',
                'format'        : 'netcdf',
                'area'          : [75.,185, 15., 320.], #N,W,S,E
                'date': dates[0].strftime("%Y-%m-%d")+"/"+dates[-1].strftime("%Y-%m-%d"),
                'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
            }
        r3 = c3.retrieve(name, request, None)

        with fsspec.open(r3.location) as f3:
            apcp_var = xr.open_dataset(f3, engine='scipy')
        #    print(apcp_var)
        
    apcp_var.to_netcdf("apcp_var_"+case_date.strftime("%Y%m%d")+".nc")

    ## get DT theta - but for this we need to go to the ECMWF queue which may take a while...

    if plot_dt_theta:

        ### create string with the times we need
        time_list = dates.strftime("%H:%M").unique().tolist()
        delim = '/'

        print('getting DT data from archive')

        c4 = cdsapi.Client()
        name = 'reanalysis-era5-complete'
        request = {
            'class': 'ea',
            'date': dates[0].strftime("%Y-%m-%d")+"/to/"+dates[-1].strftime("%Y-%m-%d"),
            'expver': '1',
            'levelist': '2000',
            'levtype': 'pv',
            'param': '3.128/54.128/131.128/132.128',
            'stream': 'oper',
            'grid'  : '0.25/0.25',
            'area'  : '75/185/15/320', #N,W,S,E
            'time': delim.join(time_list),             
            'type': 'an',
            'format'        : 'netcdf'
            }
        r4 = c4.retrieve(name, request, None)

        with fsspec.open(r4.location) as f4:
            dt_vars = xr.open_dataset(f4, engine='scipy')
        #    print(dt_vars)
        
        ## save the dt theta file for this case for use later
        dt_vars.to_netcdf("dt_vars_"+case_date.strftime("%Y%m%d")+".nc")

    print('done with downloads, moving on')
    

    #### now go on to make the maps!!

    for this_time in pd.to_datetime(iso_vars.time.values):

        iso_vars_this = iso_vars.sel(time=this_time)
        sl_vars_this = sl_vars.sel(time=this_time)
        vtime = pd.to_datetime(iso_vars_this.time.values)

        print("working on "+vtime.strftime("%H%M UTC %d %b %Y"))

        ## pull out needed variables
        lon = iso_vars_this.longitude.data
        lat = iso_vars_this.latitude.data

        if plot_500vort:

            ## get 500-mb data to start
            hght_500 = iso_vars_this.z.sel(level=500)/9.81  ## this was geopotential, convert to geopotential height
            uwnd_500 = iso_vars_this.u.sel(level=500)
            vwnd_500 = iso_vars_this.v.sel(level=500)

            # let's smooth the height a bit
            hght_500 = ndimage.gaussian_filter(hght_500, sigma=1, order=0) * units('m')

            ### calculate absolute vorticity
            avor = mpcalc.absolute_vorticity(uwnd_500,vwnd_500)

            #smooth
            avor = ndimage.gaussian_filter(avor, sigma=3, order=0) * units('1/s')
            avor = avor.squeeze()  ## remove first empty dimension

            ## vorticity advection
            #vort_adv = mpcalc.advection(avor, [uwnd_500, vwnd_500], (dx,dy)) * 1e9

            ## give units to winds & smooth slightly
            uwnd_500 = ndimage.gaussian_filter(uwnd_500, sigma=1, order=0) * units('m/s')
            vwnd_500 = ndimage.gaussian_filter(vwnd_500, sigma=1, order=0) * units('m/s')

            #Map Creation

            print('plotting 500 vorticity')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 500-hPa heights, absolute vorticity, winds',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([217.5, 290., 18., 66.])
            #ax.set_extent([227.5, 290., 20., 58.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Height Contours
            clev500 = np.arange(5100, 6061, 60)
            cs = ax.contour(lon, lat, hght_500, clev500, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                       rightside_up=True, use_clabeltext=True)

            # Absolute Vorticity colors
            # Use two different colormaps from matplotlib and combine into one color set
            clevs_500_avor = list(range(-8, 1, 1))+list(range(12, 36, 1))
            colors1 = plt.cm.YlOrRd(np.linspace(0, 1, 36))
            colors2 = plt.cm.BuPu(np.linspace(0.5, 0.75, 8))
            colors = np.vstack((colors2, (1, 1, 1, 1), colors1))

            # Plot absolute vorticity values (multiplying by 10^5 to scale appropriately)
            cf = ax.contourf(lon, lat, avor*1e5, clevs_500_avor, colors=colors, extend='both',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('absolute vorticity ($*10^5 s^{-1}$)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 14)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_500.to(units('kt'))[wind_slice,wind_slice].m, vwnd_500.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/500vort_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)


        if plot_850temp:

            ## now go on to 850-hPa temperature
            ## with some guidance from here: https://unidata.github.io/python-training/gallery/850hpa_temperature_advection/

            hght_850 = iso_vars_this.z.sel(level=850)/9.81  ## this was geopotential, convert to geopotential height
            uwnd_850 = iso_vars_this.u.sel(level=850)
            vwnd_850 = iso_vars_this.v.sel(level=850)
            tmpk_850 = iso_vars_this.t.sel(level=850)

            ## smooth the temperature 
            tmpk_850 = ndimage.gaussian_filter(tmpk_850, sigma=1, order=0) * units('kelvin')
            tmpc_850 = tmpk_850.to(units('degC'))
            hght_850 = ndimage.gaussian_filter(hght_850, sigma=1, order=0) * units('m')

            ## give units to winds & smooth slightly
            uwnd_850 = ndimage.gaussian_filter(uwnd_850, sigma=1, order=0) * units('m/s')
            vwnd_850 = ndimage.gaussian_filter(vwnd_850, sigma=1, order=0) * units('m/s')

            #Map Creation

            print('plotting 850 temperature')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 850-hPa heights, temperature, winds',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Height Contours
            clev850 = np.arange(0, 3001, 30)
            cs = ax.contour(lon, lat, hght_850, clev850, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                       rightside_up=True, use_clabeltext=True)

            # Plot colorfill and dashed contours of 850-hPa temperatures in Celsius
            clevs_850_tmpc = np.arange(-36, 38, 2)
            cf = ax.contourf(lon, lat, tmpc_850, clevs_850_tmpc, cmap=plt.cm.seismic, extend='both',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('Temperature (C)', fontsize=13)

            ## plot the 0C line
            clev_thk = [0]
            cs4 = ax.contour(lon, lat, tmpc_850, clev_thk, colors='blue',
                             linewidths=2, transform=ccrs.PlateCarree())
            plt.clabel(cs4, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 10)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_850.to(units('kt'))[wind_slice,wind_slice].m, vwnd_850.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/850temp_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        ## and now q-vectors
        if plot_700qvect:

            hght_700 = iso_vars_this.z.sel(level=700)/9.81  ## this was geopotential, convert to geopotential height
            uwnd_700 = iso_vars_this.u.sel(level=700)
            vwnd_700 = iso_vars_this.v.sel(level=700)
            tmpk_700 = iso_vars_this.t.sel(level=700)

            ## first, smooth everything pretty heavily:
            #tmpk_700 = ndimage.gaussian_filter(tmpk_700, sigma=6, order=0) * units('kelvin')
            #uwnd_700 = ndimage.gaussian_filter(uwnd_700, sigma=6, order=0) * units('m/s')
            #vwnd_700 = ndimage.gaussian_filter(vwnd_700, sigma=6, order=0) * units('m/s') 
            tmpk_700 = mpcalc.smooth_gaussian(tmpk_700,35)
            uwnd_700 = mpcalc.smooth_gaussian(uwnd_700,35) 
            vwnd_700 = mpcalc.smooth_gaussian(vwnd_700,35) 

            # less smoothing for the height
            hght_700 = ndimage.gaussian_filter(hght_700, sigma=2, order=0) * units('m')

            # Compute the Q-vector components
            uqvect, vqvect = mpcalc.q_vector(uwnd_700, vwnd_700, tmpk_700, 700*units.hPa)

            # Compute the divergence of the Q-vectors calculated above
            q_div = -2*mpcalc.divergence(uqvect, vqvect)

            ## smooth it too
            q_div = ndimage.gaussian_filter(q_div, sigma=4, order=0) * units('meter/(kilogram second)')

            #Map Creation

            print('plotting 700 q-vectors')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 700-hPa heights, temp, Q-vectors and Q-vector Conv.',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Height Contours
            clev700 = np.arange(2100, 18001, 30)
            cs = ax.contour(lon, lat, hght_700, clev700, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot temperature Contours
            clevtemp700 = np.arange(-30, 20, 2)
            cs2 = ax.contour(lon, lat, tmpk_700.metpy.convert_units('degC'), clevtemp700, colors='grey',
                             linewidths=1.5, linestyles='dashed', transform=ccrs.PlateCarree())
            plt.clabel(cs2, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)

            # Plot 700-hPa Q-Vector Divergence and scale
            clevs_qdiv = list(range(-30, -4, 5))+list(range(5, 31, 5))
            cf = ax.contourf(lon, lat, q_div*1e18, clevs_qdiv, extend='both',
                             cmap='bwr', transform=ccrs.PlateCarree())
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9,
                              ticks=clevs_qdiv)
            cb.set_label('Q-Vector Conv. (2*10$^{18}$ m s$^{-1}$ kg$^{-1}$)', fontsize=13)

            # Plot 700-hPa Q-vectors, scale to get nice sized arrows
            wind_slice = slice(None, None, 10)
            ax.quiver(lon[wind_slice],lat[wind_slice],uqvect[wind_slice,wind_slice].values,vqvect[wind_slice,wind_slice].values,
                      pivot='mid', color='black',
                      scale=1e-11, scale_units='inches',
                      transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/700qvect_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        ## do 700-mb wind speed also

        if plot_700wsp:

            hght_700 = iso_vars_this.z.sel(level=700)/9.81  ## this was geopotential, convert to geopotential height
            uwnd_700 = iso_vars_this.u.sel(level=700)
            vwnd_700 = iso_vars_this.v.sel(level=700)
            tmpk_700 = iso_vars_this.t.sel(level=700)

            # Use MetPy to calculate the wind speed for colorfill plot
            sped_700 = mpcalc.wind_speed(uwnd_700, vwnd_700)
            sped_700 = ndimage.gaussian_filter(sped_700, sigma=2, order=0) * units('m/s')

            ## smooth the height, u and v
            hght_700 = ndimage.gaussian_filter(hght_700, sigma=1, order=0) * units('m')
            uwnd_700 = ndimage.gaussian_filter(uwnd_700, sigma=1, order=0) * units('m/s')
            vwnd_700 = ndimage.gaussian_filter(vwnd_700, sigma=1, order=0) * units('m/s')

            #Map Creation

            print('plotting 700 wind speed')

            fig = plt.figure(1, figsize=(20., 16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 700-hPa heights, winds, and isotachs',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            #ax.set_extent([232.5, 285., 20., 52.5])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Height Contours
            clev700 = np.arange(2100, 18001, 30)
            cs = ax.contour(lon, lat, hght_700, clev700, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot Colorfill Wind Speeds in m/s
            clevs_700_sped = np.arange(15, 55, 5)
            cf = ax.contourf(lon, lat, sped_700, clevs_700_sped, cmap=plt.cm.YlOrRd,extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('wind speed (m/s)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 10)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_700.to(units('kt'))[wind_slice,wind_slice].m, vwnd_700.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/700wsp_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)


        ## now go on to 250-hPa isotachs
        if plot_250wsp:

            hght_250 = iso_vars_this.z.sel(level=250)/9.81  ## this was geopotential, convert to geopotential height
            uwnd_250 = iso_vars_this.u.sel(level=250)
            vwnd_250 = iso_vars_this.v.sel(level=250)
            tmpk_250 = iso_vars_this.t.sel(level=250)

            # Use MetPy to calculate the wind speed for colorfill plot
            sped_250 = mpcalc.wind_speed(uwnd_250, vwnd_250)
            sped_250 = ndimage.gaussian_filter(sped_250, sigma=2, order=0) * units('m/s')

            ## smooth the height, u and v
            hght_250 = ndimage.gaussian_filter(hght_250, sigma=1, order=0) * units('m')
            uwnd_250 = ndimage.gaussian_filter(uwnd_250, sigma=1, order=0) * units('m/s')
            vwnd_250 = ndimage.gaussian_filter(vwnd_250, sigma=1, order=0) * units('m/s')

            #Map Creation

            print('plotting 250 wind speed')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 250-hPa heights, winds, and isotachs',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([217.5, 290., 18., 66.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Height Contours
            clev250 = np.arange(6000, 21001, 120)
            cs = ax.contour(lon, lat, hght_250, clev250, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot Colorfill Wind Speeds in m/s
            clevs_250_sped = np.arange(30, 110, 10)
            cf = ax.contourf(lon, lat, sped_250, clevs_250_sped, cmap=plt.cm.YlOrRd,extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('wind speed (m/s)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 14)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_250.to(units('kt'))[wind_slice,wind_slice].m, vwnd_250.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/250wsp_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        ## Now MSLP and 700-mb RH
        if plot_mslp_thk:

            rh_700 = iso_vars_this.r.sel(level=700)

            hght_1000 = iso_vars_this.z.sel(level=1000)/9.81
            hght_1000 = ndimage.gaussian_filter(hght_1000, sigma=1, order=0) * units('m')
            hght_500 = iso_vars_this.z.sel(level=500)/9.81
            hght_500 = ndimage.gaussian_filter(hght_500, sigma=1, order=0) * units('m')

            thk_1000_500 = hght_500 - hght_1000

            mslp = ndimage.gaussian_filter(sl_vars_this.msl.data, sigma=2, order=0) * units('Pa')
            thk_1000_500 = ndimage.gaussian_filter(thk_1000_500, sigma=1, order=0) * units('m')

            #Map Creation

            print('plotting MSLP/RH')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title("ERA5 MSLP, 1000{0}500".format(enDash)+"-hPa thickness, and 700-hPa RH",
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot MSLP Contours
            clev_mslp = np.arange(940,1060,4)
            cs = ax.contour(lon, lat, mslp.to(units('hPa')), clev_mslp, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot thickness Contours
            clev_thk = np.arange(3600,5400,60)
            cs2 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='blue',
                             linewidths=1.5, linestyles='dashed', transform=ccrs.PlateCarree())
            plt.clabel(cs2, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)
            ## above 5400
            clev_thk = np.arange(5460,6600,60)
            cs3 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='red',
                             linewidths=1.5, linestyles='dashed', transform=ccrs.PlateCarree())
            plt.clabel(cs3, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)
            ## and the 5400 line
            clev_thk = [5400]
            cs4 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='red',
                             linewidths=2, transform=ccrs.PlateCarree())
            plt.clabel(cs4, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)

            # Plot Colorfill RH
            clevs_700_rh = np.arange(70,105,5)
            cf = ax.contourf(lon, lat, rh_700, clevs_700_rh, cmap=plt.cm.Greens, extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('RH (%)', fontsize=13)

            outfile=case_date.strftime("%Y%m%d")+"/mslp_thk_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        if plot_pwat:

            print('plotting PW')

            uwnd_700 = iso_vars_this.u.sel(level=700)
            vwnd_700 = iso_vars_this.v.sel(level=700)

            ## smooth the u and v
            uwnd_700 = ndimage.gaussian_filter(uwnd_700, sigma=1, order=0) * units('m/s')
            vwnd_700 = ndimage.gaussian_filter(vwnd_700, sigma=1, order=0) * units('m/s')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 precipitable water and 700-mb winds',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Colorfill PW
            clevs_pwat = np.arange(5,70,5)
            cf = ax.contourf(lon,lat, sl_vars_this.tcwv.data, clevs_pwat, cmap=plt.cm.Greens,extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('precipitable water (mm)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 10)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_700.to(units('kt'))[wind_slice,wind_slice].m, vwnd_700.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/pwat_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        # now precip/MSLP

        if plot_mslp_precip:

            # convert from m to mm and to inches
            apcp = 1000.*apcp_var.tp
            apcp = apcp/25.4

            ## and calculate sum over 12-hour period
            apcp_tot = apcp.rolling(time=6).sum()

            # find this time
            apcp_tot_this = apcp_tot.sel(time=this_time)

            hght_1000 = iso_vars_this.z.sel(level=1000)/9.81
            hght_1000 = ndimage.gaussian_filter(hght_1000, sigma=1, order=0) * units('m')
            hght_500 = iso_vars_this.z.sel(level=500)/9.81
            hght_500 = ndimage.gaussian_filter(hght_500, sigma=1, order=0) * units('m')

            thk_1000_500 = hght_500 - hght_1000

            mslp = ndimage.gaussian_filter(sl_vars_this.msl.data, sigma=2, order=0) * units('Pa')
            thk_1000_500 = ndimage.gaussian_filter(thk_1000_500, sigma=1, order=0) * units('m')

            # plot MSLP/thickness/precip
            #Map Creation

            print('plotting MSLP/thickness/precip')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title("ERA5 MSLP, 1000{0}500".format(enDash)+"-hPa thickness, and precip in previous 6 hours",
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([227.5, 290., 20., 58.])
            #ax.set_extent([232.5, 285., 20., 52.5])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot MSLP Contours
            clev_mslp = np.arange(940,1060,4)
            cs = ax.contour(lon, lat, mslp.to(units('hPa')), clev_mslp, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot thickness Contours
            clev_thk = np.arange(3600,5400,60)
            cs2 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='blue',
                             linewidths=1.5, linestyles='dashed', transform=ccrs.PlateCarree())
            plt.clabel(cs2, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)
            ## above 5400
            clev_thk = np.arange(5460,6600,60)
            cs3 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='red',
                             linewidths=1.5, linestyles='dashed', transform=ccrs.PlateCarree())
            plt.clabel(cs3, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)
            ## and the 5400 line
            clev_thk = [5400]
            cs4 = ax.contour(lon, lat, thk_1000_500, clev_thk, colors='red',
                             linewidths=2, transform=ccrs.PlateCarree())
            plt.clabel(cs4, fontsize=12, inline=1, inline_spacing=10, fmt='%i',
                             rightside_up=True, use_clabeltext=True)

            # Plot Colorfill precip
            # draw filled contours.
            cmap_data = [(1.000000, 1.000000, 1.000000),
            (0.498039, 1.000000, 0.000000),
            (0.000000, 0.803922, 0.000000),
            (0.000000, 0.545098, 0.000000),
            (0.062745, 0.305882, 0.545098),
            (0.117647, 0.564706, 1.000000),
            (0.000000, 0.698039, 0.933333),
            (0.000000, 0.933333, 0.933333),
            (0.537255, 0.407843, 0.803922),
            (0.568627, 0.172549, 0.933333),
            (0.545098, 0.000000, 0.545098),
            (0.545098, 0.000000, 0.000000),
            (0.803922, 0.000000, 0.000000),
            (0.933333, 0.250980, 0.000000),
            (1.000000, 0.498039, 0.000000),
            (0.803922, 0.521569, 0.000000),
            (1.000000, 0.843137, 0.000000)]
            cmap = mcolors.ListedColormap(cmap_data, 'precip_rss')

            clevs = [0,0.01,0.05,0.1,0.25,0.5,0.75,1.00,1.5,2,2.5,3,3.5,4,5,6,8]
            norm = mcolors.BoundaryNorm(clevs, cmap.N)

            # In future MetPy
            #norm, cmap = ctables.registry.get_with_boundaries('precipitation', clevs)

            cf = ax.contourf(apcp_var.longitude.data, apcp_var.latitude.data, apcp_tot_this, 
                             clevs, cmap=cmap, norm=norm, extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9,
                             ticks=clevs)
            cb.set_label('precipitation (inches)', fontsize=13)

            outfile=case_date.strftime("%Y%m%d")+"/mslp_prcp_thk_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        ## now go on to surface dewpoint/winds/MSLP

        if plot_sfc_dwpt:

            mslp = ndimage.gaussian_filter(sl_vars_this.msl.data, sigma=2, order=0) * units('Pa')
            sfc_dwpt = sl_vars_this.d2m.data - 273.15
            sfc_dwpt = ndimage.gaussian_filter(sfc_dwpt, sigma=0, order=0) * units('degC')
            uwnd_10m = ndimage.gaussian_filter(sl_vars_this.u10.data, sigma=1, order=0) * units('m/s')
            vwnd_10m = ndimage.gaussian_filter(sl_vars_this.v10.data, sigma=1, order=0) * units('m/s')

            print('plotting MSLP/sfc dewpoint/winds')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 MSLP, 2-m dewpoint, 10-m winds',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([235, 285., 23., 52.5])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot MSLP Contours
            clev_mslp = np.arange(940,1060,4)
            cs = ax.contour(lon, lat, mslp.to(units('hPa')), clev_mslp, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot colorfill and dashed contours of sfc dewpoints in Celsius
            clevs_sfc_dwpt = np.arange(0,28,2)
            cf = ax.contourf(lon, lat, sfc_dwpt, clevs_sfc_dwpt, cmap=plt.cm.Greens, extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('dewpoint (C)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 6)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_10m.to(units('kt'))[wind_slice,wind_slice].m, vwnd_10m.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=6.5, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/sfc_dwpt_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        ## now go on to MUCAPE and sfc-500 mb shear

        if plot_mucape:

            uwnd_500 = iso_vars_this.u.sel(level=500)
            vwnd_500 = iso_vars_this.v.sel(level=500)
            uwnd_10m = sl_vars_this.u10
            vwnd_10m = sl_vars_this.v10

            ushr_sfc_500 = (uwnd_500 - uwnd_10m) * units('m/s')
            vshr_sfc_500 = (vwnd_500 - vwnd_10m) * units('m/s')
            shearmag_sfc_500 = mpcalc.wind_speed(ushr_sfc_500, vshr_sfc_500)
            shearmag_sfc_500 = ndimage.gaussian_filter(shearmag_sfc_500, sigma=2, order=0) * units('m/s')

            ## smooth the vectors
            ushr_sfc_500 = ndimage.gaussian_filter(ushr_sfc_500, sigma=1, order=0) * units('m/s')
            vshr_sfc_500 = ndimage.gaussian_filter(vshr_sfc_500, sigma=1, order=0) * units('m/s')

            mucape = ndimage.gaussian_filter(sl_vars_this.cape.data, sigma=2, order=0) * units('J/kg')

            print('plotting MUCAPE/sfc-500 mb shear')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 MUCAPE, sfc-500-hPa shear (m/s)',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([235, 285., 23., 52.5])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot shear magnitude contours
            clev_bwd = np.arange(10,90,5)
            cs = ax.contour(lon, lat, shearmag_sfc_500, clev_bwd, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot colorfill and dashed contours of 850-hPa dewpoints in Celsius
            clevs_mucape = [100,250,500,1000,1500,2000,2500,3000,3500,4000]
            cf = ax.contourf(lon, lat, mucape, clevs_mucape, cmap=plt.cm.YlOrRd, extend='max',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('MUCAPE (J/kg)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 6)
            ax.barbs(lon[wind_slice], lat[wind_slice], ushr_sfc_500.to(units('kt'))[wind_slice,wind_slice].m, vshr_sfc_500.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/mucape_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        if plot_sfc_rh:

            mslp = ndimage.gaussian_filter(sl_vars_this.msl.data, sigma=2, order=0) * units('Pa')
            sfc_dwpt = sl_vars_this.d2m
            sfc_temp = sl_vars_this.t2m
            sfc_rh = mpcalc.relative_humidity_from_dewpoint(sfc_temp, sfc_dwpt)

            sfc_rh = ndimage.gaussian_filter(sfc_rh, sigma=1, order=0) * units('degC')        
            uwnd_10m = ndimage.gaussian_filter(sl_vars_this.u10.data, sigma=1, order=0) * units('m/s')
            vwnd_10m = ndimage.gaussian_filter(sl_vars_this.v10.data, sigma=1, order=0) * units('m/s')

            print('plotting MSLP/sfc RH/winds')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 MSLP, 2-m relative humidity, 10-m winds',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([235, 285., 23., 52.5])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot MSLP Contours
            clev_mslp = np.arange(940,1060,4)
            cs = ax.contour(lon, lat, mslp.to(units('hPa')), clev_mslp, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot colorfill and dashed contours of RH
            clevs_sfc_rh = np.arange(0,1.05,0.05)
            cf = ax.contourf(lon, lat, sfc_rh, clevs_sfc_rh, cmap=BrBG_11.mpl_colormap, extend='both',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('relative humidity (%)', fontsize=13)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 6)
            ax.barbs(lon[wind_slice], lat[wind_slice], uwnd_10m.to(units('kt'))[wind_slice,wind_slice].m, vwnd_10m.to(units('kt'))[wind_slice,wind_slice].m, 
                     length=6.5, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/sfc_rh_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)

        # Now DT theta and winds

        if plot_dt_theta:

            dt_vars_this = dt_vars.sel(time=this_time)

            ## smooth the theta and winds
            dt_theta = ndimage.gaussian_filter(dt_vars_this.pt.data, sigma=3, order=0) * units('kelvin')
            dt_uwnd = ndimage.gaussian_filter(dt_vars_this.u.data, sigma=2, order=0) * units('kt')
            dt_vwnd = ndimage.gaussian_filter(dt_vars_this.v.data, sigma=2, order=0) * units('kt')

            dt_lon = dt_vars_this.longitude.data
            dt_lat = dt_vars_this.latitude.data

            ## also calculate the low-level relative vort:
            uwnd_850 = iso_vars_this.u.sel(level=850)
            vwnd_850 = iso_vars_this.v.sel(level=850)
            relvor_850 = mpcalc.vorticity(uwnd_850,vwnd_850) 

            #smooth
            relvor_850 = ndimage.gaussian_filter(relvor_850, sigma=3, order=0) * units('1/s')
            relvor_850 = relvor_850.squeeze() ## remove first empty dimension

            #Map Creation

            print('plotting DT theta')

            fig = plt.figure(1, figsize=(20,16))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                                   hspace=0.01, wspace=0.01)
            ax = plt.subplot(gs[0], projection=plotcrs)

            # Plot Titles
            plt.title(r'ERA5 dynamic tropopause potential temperature and winds; 850-hPa relative vorticity',
                      loc='left', fontsize=16)
            plt.title(vtime.strftime("%H%M UTC %d %b %Y"), loc='right', fontsize=16)

            # Plot Background
            ax.set_extent([217.5, 290., 18., 66.])
            ax.coastlines('50m', edgecolor='black', linewidth=0.75)
            ax.add_feature(cfeature.STATES, linewidth=.5)

            # Plot Colorfill DT theta
            clevs_dt = np.arange(260,400,5)
            cf = ax.contourf(dt_lon, dt_lat, dt_theta, clevs_dt, cmap=plt.cm.twilight_shifted,extend='both',
                             transform=datacrs)
            cb = plt.colorbar(cf, orientation='horizontal', pad=0.01, aspect=50, shrink=0.9)
            cb.set_label('potential temperature (K)', fontsize=13)

            # Plot 850-mb vort contours
            clev_relvor = np.arange(4,40,4)
            cs = ax.contour(lon, lat, relvor_850*1e5, clev_relvor, colors='black', linewidths=2,
                            linestyles='solid', transform=ccrs.PlateCarree())
            plt.clabel(cs, fontsize=12, inline=1, inline_spacing=10, fmt='%i', 
                            rightside_up=True, use_clabeltext=True)

            # Plot Wind Barbs
            # Transform Vectors and plot wind barbs.
            wind_slice = slice(None, None, 14)
            ax.barbs(dt_lon[wind_slice], dt_lat[wind_slice], dt_uwnd[wind_slice,wind_slice].m, dt_vwnd[wind_slice,wind_slice].m, 
                     length=7, linewidth=0.5,pivot='middle', transform=ccrs.PlateCarree())

            outfile=case_date.strftime("%Y%m%d")+"/dt_theta_era5_"+this_time.strftime('%Y%m%d')+"_"+this_time.strftime('%H%M')+".png"
            fig.savefig(outfile)
            plt.close(fig)



# In[ ]:





# In[ ]:




