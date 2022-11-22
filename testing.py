#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:55:47 2022

@author: thoverga
"""


#example see: https://developers.google.com/earth-engine/tutorials/community/intro-to-python-api


import ee
import pandas as pd
import datetime
import folium
import pprint
#%%

ee.Authenticate()
ee.Initialize()
#%%

# Import the MODIS land cover collection.
lc = ee.ImageCollection('MODIS/006/MCD12Q1')
# See link: https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD12Q1

# Import the MODIS land surface temperature collection.
lst = ee.ImageCollection('MODIS/006/MOD11A1')

# Import the USGS ground elevation image.
elv = ee.Image('USGS/SRTMGL1_003')

# All of these images come in a different resolution, frequency, and possibly projection, ranging from daily images in a 1 km resolution for LST (hence an ee.ImageCollection — a collection of several ee.Images) to a single image representing data for the year 2000 in a 30 m resolution for the ELV. While we need to have an eye on the frequency, GEE takes care of resolution and projection by resampling and reprojecting all data we are going to work with to a common projection (learn more about projections in Earth Engine). We can define the resolution (called scale in GEE) whenever necessary and of course have the option to force no reprojection.

# As you can see in the description of the datasets, they include several sets of information stored in several bands. For example, these bands are associated with the LST collection:

# LST_Day_1km: Daytime Land Surface Temperature
# Day_view_time: Local time of day observation
# LST_Night_1km: Nighttime Land Surface Temperature
# etc.
# The description page of the collection tells us that the name of the band associated with the daytime LST is LST_Day_1km which is in units of Kelvin. In addition, values are ranging from 7,500 to 65,535 with a corrective scale of 0.02.

# Then, we have to filter the collection on the period of time we want. We can do that using the filterDate() method. We also need to select the bands we want to work with. Therefore, we decide to focus on daytime LST so we select the daytime band LST_Day_1km and its associated quality indicator QC_Day with the select() method.

# Initial date of interest (inclusive).
i_date = '2017-01-01'

# Final date of interest (exclusive).
f_date = '2020-01-01'

# Selection of appropriate bands and dates for LST.
lst = lst.select('LST_Day_1km', 'QC_Day').filterDate(i_date, f_date)

# Define the urban location of interest as a point near Lyon, France.
u_lon = 4.8148
u_lat = 45.7758
u_poi = ee.Geometry.Point(u_lon, u_lat)

# Define the rural location of interest as a point away from the city.
r_lon = 5.175964
r_lat = 45.574064
r_poi = ee.Geometry.Point(r_lon, r_lat)

# =============================================================================
# Scale 
# =============================================================================
scale = 1000  # scale in meters


# # =============================================================================
# # Point extractor
# # =============================================================================

# # Print the elevation near Lyon, France.
# elv_urban_point = elv.sample(u_poi, scale).first().get('elevation').getInfo()
# print('Ground elevation at urban point:', elv_urban_point, 'm')


# # =============================================================================
# # #Aggregation functins
# # =============================================================================
# # Calculate and print the mean value of the LST collection at the point.
# lst_urban_point = lst.mean().sample(u_poi, scale).first().get('LST_Day_1km').getInfo()
# print('Average daytime LST at urban point:', round(lst_urban_point*0.02 -273.15, 2), '°C')

# # Print the land cover type at the point.
# lc_urban_point = lc.first().sample(u_poi, scale).first().get('LC_Type1').getInfo()
# print('Land cover value at urban point is:', lc_urban_point)



# =============================================================================
# Data extraction
# =============================================================================

def extract_pointvalue(dataset, lat, lon, band=None, scale=1000):
    
    raster = ee.Image(dataset)
    
    if isinstance(band, type(None)):
        bandnames = raster.bandNames().getInfo()
        
        assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
        band = bandnames[0]
        
    
    
    location = ee.Geometry.Point(lon, lat)
    
    data = raster.sample(location, scale).first().get(band).getInfo()
    
    return data



landcover_at_point = extract_pointvalue(dataset=lc,
                                        lat = u_lat,
                                        lon = u_lon)

print(landcover_at_point)





def extract_timeseries_at_point(dataset, band, lat, lon, startdatestr, enddatestr, scale=1000):
    
    dataset = ee.ImageCollection(dataset).select(band).filterDate(ee.Date(i_date), ee.Date(f_date))
    
    location = ee.Geometry.Point(lon, lat)
    
    data = dataset.getRegion(location, scale).getInfo()
    
    df = pd.DataFrame(data=data[1:], columns=data[0])
    
    return df

# test = extract_timeseries_at_point(dataset = "ECMWF/ERA5_LAND/MONTHLY_BY_HOUR",
#                                    band = 'temperature_2m',
#                                    lat = 51.1662 ,
#                                    lon = 4.3665,
#                                    startdatestr = '2019-01-01',
#                                    enddatestr = '2020-01-01',
#                                    scale=1000)

# test['temperature_2m'].plot()





def get_info(obj):
    obj_info = dict(ee.Algorithms.Describe(obj).getInfo())
    print('--------------  INFO  ----------------')
    print(obj_info['properties']['description'])    
    print('-------------- data INFO  ----------------')
    print(obj_info['id'])
    print(obj_info['bands'])

get_info(elv)



#%% plot region box
import folium
import branca.colormap as cm


def make_ee_plot(ee_datasetname,band=None, left_bottom_lon=2.4, left_bottom_lat=50.8,
                 right_top_lon=5.61, righ_top_lat=51.45,
                 vmin=None, vmax=None, mincolor='white', maxcolor='blue',
                 outputfile="/home/thoverga/Documents/github/google_earth_engine/mymap.html"):

    
# =============================================================================
#     helper function
# =============================================================================
    def add_ee_layer(self, ee_image_object, vis_params, name):
      map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
      folium.raster_layers.TileLayer(
          tiles=map_id_dict['tile_fetcher'].url_format,
          attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
          name=name,
          overlay=True,
          control=True
      ).add_to(self)
    
    folium.Map.add_ee_layer = add_ee_layer

# =============================================================================
# create Roi rectangle
# =============================================================================

    
    pointlist = [ee.Geometry.Point(left_botom_lon, left_botom_lat), 
                 ee.Geometry.Point(right_top_lon, right_top_lat)]
    
    rectangle = ee.Geometry.Rectangle(pointlist)

# =============================================================================
# get gee data and clip to roi
# =============================================================================
   
    dataset = ee.Image(ee_datasetname)
    
    if isinstance(band, type(None)):
        bandnames = dataset.bandNames().getInfo()
        
        assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
        band = bandnames[0]
    
     
    clipped_raster = dataset.clip(rectangle)

   

# =============================================================================
# vizualisation parameters
# =============================================================================
    if isinstance(vmin, type(None)):
        vmin=clipped_raster.reduceRegion(ee.Reducer.min(), rectangle).getInfo()[band]
    if isinstance(vmax, type(None)):
        vmax = clipped_raster.reduceRegion(ee.Reducer.max(), rectangle).getInfo()[band]
    image_viz_params = {
        'bands': band,
        'min': vmin,
        'max': vmax,
        'palette': [mincolor, maxcolor]
        }
    
    # Define a colormap.
    rech_colormap = cm.LinearColormap(
        colors=image_viz_params["palette"],
        vmin=image_viz_params["min"],
        vmax=image_viz_params["max"],
        )


# =============================================================================
# Build map layers
# =============================================================================

    
    mapfig = folium.Map()
    mapfig.add_ee_layer(clipped_raster, image_viz_params, 'testnaam')
    
    mapfig.fit_bounds(bounds=[(left_botom_lat, left_botom_lon), (right_top_lat, right_top_lon)])
    
    mapfig.add_child(rech_colormap)

# =============================================================================
# save map
# =============================================================================
    mapfig.save(outputfile)


make_ee_plot(ee_datasetname = 'USGS/SRTMGL1_003', vmin=-30, vmax=60)



#%%
import geemap
Map = geemap.Map()


dem_vis = {
    'bands': band,
    'min': clipped_raster.reduceRegion(ee.Reducer.min(), rectangle).getInfo()[band],
    'max': clipped_raster.reduceRegion(ee.Reducer.max(), rectangle).getInfo()[band],
}


Map.addLayer(clipped_raster, dem_vis)

Map.save("/home/thoverga/Documents/github/google_earth_engine/mymap.html")

#%% Helper functions
def datetime_to_ee_format(dt):
    return ee.Date(dt, opt_tz='UTC')
    
def ee_to_datetime_format(ee_date):
   
    py_date = datetime.datetime.utcfromtimestamp(ee_date.getInfo()['value']/1000.0)
    print(py_date)

start_period =  datetime.datetime(2020, 7, 1)
test = datetime_to_ee_form(start_period)

# ee_to_datetime_format(test)


#%% Era5 extractor

dataset_name = "ECMWF/ERA5_LAND/HOURLY"

start_period =  datetime.datetime(2020, 7, 1)
end_period = datetime.datetime(2020, 7, 1, 10)



loc_lat = 51.174
loc_lon = 4.104



#convert to ee classes
ee_startdt = datetime_to_ee_format(start_period)
ee_enddt = datetime_to_ee_format(end_period)
ee_locloc = ee.Geometry.Point([loc_lat, loc_lon])


#extract from dataset
# image = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY").filter(ee.Filter.date('2020-07-01', '2020-07-02'))

# visualization = {
#   'bands': ['temperature_2m'],
#   'min': 250.0,
#   'max': 320.0,
#   'palette': [
#     "#000080","#0000D9","#4000FF","#8000FF","#0080FF","#00FFFF",
#     "#00FF80","#80FF00","#DAFF00","#FFFF00","#FFF500","#FFDA00",
#     "#FFB000","#FFA400","#FF4F00","#FF2500","#FF0A00","#FF00FF",
#   ]
# }
#%%
import folium


def add_ee_layer(self, ee_image_object, vis_params, name):
  map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
  folium.raster_layers.TileLayer(
      tiles=map_id_dict['tile_fetcher'].url_format,
      attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
      name=name,
      overlay=True,
      control=True
  ).add_to(self)

folium.Map.add_ee_layer = add_ee_layer

# Load an image.
image = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_044034_20140318')

# Define the visualization parameters.
image_viz_params = {
    'bands': ['B5', 'B4', 'B3'],
    'min': 0,
    'max': 0.5,
    'gamma': [0.95, 1.1, 1]
}

# Define a map centered on San Francisco Bay.
map_l8 = folium.Map(location=[37.5010, -122.1899], zoom_start=10)

# Add the image layer to the map and display it.
map_l8.add_ee_layer(image, image_viz_params, 'false color composite')
display(map_l8)

map_l8.save("/home/thoverga/Documents/github/google_earth_engine/mymap.html")


#%% handy functions
def show_collection_date_range(imagecol):
    dtrange = imagecol.reduceColumns(ee.Reducer.minMax(), ['system:time_start'])
    print('Date range: ', dtrange)


test = show_collection_date_range(ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")).


#%%

test = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY").filterDate(ee_startdt, ee_enddt).filterBounds(ee_locloc)