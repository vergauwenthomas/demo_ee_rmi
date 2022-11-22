#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 09:41:16 2022

@author: thoverga
"""

import pandas as pd
import geopandas as gpd
from datetime import datetime
import gee_functions as gee






#%%
# =============================================================================
# Make shure to connect to the google servers 
# =============================================================================
gee.connect_to_gee()







#%%
# =============================================================================
# Demo 1. DEM and the RMI
# =============================================================================

#extract point value

Nasa_DEM = "USGS/SRTMGL1_003" # https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003



gee.get_info(Nasa_DEM)


#RMI
RMI_lon = 4.358340
RMI_lat = 50.797343

scale = 500 #get map at 500m resolution 
height_RMI = gee.extract_pointvalue(datasetname=Nasa_DEM,
                              lat=RMI_lat,
                              lon=RMI_lon,
                              scale=scale)

print (f'The height of the RMI is {height_RMI}m, for the DEM at spatial resultion of {scale}m')


#Multiple points:
    
height_points = gee.extract_multiple_pointvalues(datasetname=Nasa_DEM,
                                                 lats=[45.7758, 45.7858],
                                                 lons=[4.8418, 4.8418])

print(f'Height of points: {height_points}')




















#%%
# =============================================================================
# Demo 2. Timeseries at the RMI
# =============================================================================


#Dataset
ERA5 = "ECMWF/ERA5/DAILY"


bandnames = ["mean_2m_air_temperature", "u_component_of_wind_10m"]

#RMI
RMI_lon = 4.358340
RMI_lat = 50.797343

#period
Startperiod='1995-01-01'
Endperiod = '2005-01-01'

timeseries_RMI = gee.extract_timeseries_at_point(datasetname=ERA5,
                                             lat=RMI_lat,
                                             lon=RMI_lon,
                                             startdatestr=Startperiod,
                                             enddatestr=Endperiod,
                                             list_of_bands=bandnames,
                                             scale=2500)


print(timeseries_RMI.head())



timeseries_RMI['mean_2m_air_temperature'].plot()












#%%
# =============================================================================
# Demo 3. Spatial Plot 
# =============================================================================
import ee


# 1. Spatial plot of IMAGE (static)


raster = "USGS/SRTMGL1_003" #DEM NASA


gee.make_ee_plot(ee_dataset=raster,
                 band=None, #If one band available, use that band
                 left_bottom_lon=2.4, #Flanders
                 left_bottom_lat=50.8, #Flanders
                 right_top_lon=5.61, #Flanders
                 right_top_lat=51.45, #Flanders
                 vmin=0,
                 vmax=150,
                 mincolor='white',
                 maxcolor='blue',
                 outputfile="/home/thoverga/Documents/github/demo_ee_rmi/mymap.html"
                 )









#%%
# 2. Spatial plot of IMAGECOLLECTION (time depending)

ERA5 = "ECMWF/ERA5/DAILY"

band = "mean_2m_air_temperature"

moment = '2010-10-20'

#convert ImageCollection to image for the given band and timeinstance
raster = ee.ImageCollection(ERA5).select(band).filterDate(moment).first()

gee.make_ee_plot(ee_dataset=raster,
                 band=band, 
                 left_bottom_lat=50.8, #Flanders
                 right_top_lon=5.61, #Flanders
                 right_top_lat=51.45, #Flanders
                 vmin=277,
                 vmax=282,
                 mincolor='blue',
                 maxcolor='red',
                 outputfile="/home/thoverga/Documents/github/demo_ee_rmi/mymap.html"
                 )


