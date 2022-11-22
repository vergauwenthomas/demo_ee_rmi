#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 09:30:15 2022

@author: thoverga
"""

import ee
import pandas as pd
import folium
import branca.colormap as cm

#TODO: documentation
#TODO: image, imagecollection try-exept statements
#TODO: imagecollection to image filters (xarray.sel()-like)

# =============================================================================
#  Connection functions
# =============================================================================

def connect_to_gee():
    ee.Authenticate()
    ee.Initialize()
    return


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


def ee_array_to_df(arr, list_of_bands):
    """Transforms client-side ee.Image.getRegion array to pandas.DataFrame."""
    df = pd.DataFrame(arr)

    # Rearrange the header.
    headers = df.iloc[0]
    df = pd.DataFrame(df.values[1:], columns=headers)

    # Remove rows without data inside.
    df = df[['longitude', 'latitude', 'time', *list_of_bands]].dropna()

    # Convert the data to numeric values.
    for band in list_of_bands:
        df[band] = pd.to_numeric(df[band], errors='coerce')

    # Convert the time field into a datetime.
    df['datetime'] = pd.to_datetime(df['time'], unit='ms')

    # Keep the columns of interest.
    df = df[['time','datetime',  *list_of_bands]]

    return df
# =============================================================================
#   Get info
# =============================================================================

def get_info(obj):
    if isinstance(obj, str):
        try: 
            obj_info = ee.Algorithms.Describe(ee.Image(obj)).getInfo()
        except ee.ee_exception.EEException:
            try:
                obj_info = ee.Algorithms.Describe(ee.ImageCollection(obj)).getInfo()
            except ee.ee_exception.EEException:
                print(obj, ' is not found as image or imagecollection!!')
                return
    else:
        obj_info = dict(ee.Algorithms.Describe(obj).getInfo())
    print('--------------  INFO  ----------------')
    print(obj_info['properties']['description'])    
    print('-------------- data INFO  ----------------')
    print(obj_info['id'])
    print(obj_info['bands'])


# =============================================================================
#   Data extraction 
# =============================================================================
def extract_pointvalue(datasetname, lat, lon, band=None, scale=1000):
    
    raster = ee.Image(datasetname)
    
    if isinstance(band, type(None)):
        bandnames = raster.bandNames().getInfo()
        
        assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
        band = bandnames[0]
        
    
    
    location = ee.Geometry.Point(lon, lat)
    
    data = raster.sample(location, scale).first().get(band).getInfo()
    
    return data

def extract_multiple_pointvalues(datasetname, lats, lons, band=None, scale=1000):
        
    features_list = []
    for point in list(zip(lons, lats)):
        features_list.append(ee.Feature(ee.Geometry.Point(point[0], point[1])))

    
    raster = ee.Image(datasetname)
    
    if isinstance(band, type(None)):
        bandnames = raster.bandNames().getInfo()
        
        assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
        band = bandnames[0]
        
    
    data = raster.sampleRegions(collection=features_list,scale=scale).getInfo()
    

    return [feat['properties'][band] for feat in data['features']]


def extract_timeseries_at_point(datasetname, lat, lon, startdatestr, enddatestr,list_of_bands=None, scale=1000):
    
    # raster = ee.Image(datasetname)
    # if isinstance(list_of_bands, type(None)):
    #     bandnames = raster.bandNames().getInfo()
        
    #     assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
    #     band = bandnames[0]
    
    
    dataset = ee.ImageCollection(datasetname).select(list_of_bands).filterDate(ee.Date(startdatestr), ee.Date(enddatestr))
    
    location = ee.Geometry.Point(lon, lat)
    
    data = dataset.getRegion(location, scale).getInfo()
    
    #array to dataframe
    df = ee_array_to_df(arr=data, list_of_bands=list_of_bands)
    
    return df

# =============================================================================
#   Spatial  plotting
# =============================================================================

def make_ee_plot(ee_dataset,
                 band=None, left_bottom_lon=2.4, left_bottom_lat=50.8,
                 right_top_lon=5.61, right_top_lat=51.45,
                 vmin=None, vmax=None, mincolor='white', maxcolor='blue',
                 outputfile="/home/thoverga/Documents/github/google_earth_engine/mymap.html"):

    


# =============================================================================
# create Roi rectangle
# =============================================================================

    
    pointlist = [ee.Geometry.Point(left_bottom_lon, left_bottom_lat), 
                 ee.Geometry.Point(right_top_lon, right_top_lat)]
    
    rectangle = ee.Geometry.Rectangle(pointlist)

# =============================================================================
# get gee data and clip to roi
# =============================================================================
    if isinstance(ee_dataset, ee.image.Image):
        dataset = ee_dataset
    else:
        dataset = ee.Image(ee_dataset)
    
    if isinstance(band, type(None)):
        bandnames = dataset.bandNames().getInfo()
        
        assert len(bandnames) == 1, 'Multiple bands are availible, specify one.'
        band = bandnames[0]
    
     
    clipped_raster = dataset.clip(rectangle)

   

# =============================================================================
# vizualisation parameters
# =============================================================================
    if isinstance(vmin, type(None)):
        # print(clipped_raster.reduceRegion(ee.Reducer.min(), rectangle).getInfo())
        # vmin=clipped_raster.reduceRegion(ee.Reducer.min(), rectangle).getInfo()[band]
        vmin=clipped_raster.reduceRegion(reducer = ee.Reducer.min(),
                                         geometry = rectangle,
                                         bestEfford=True).getInfo()[band]
    if isinstance(vmax, type(None)):
        # vmax = clipped_raster.reduceRegion(ee.Reducer.max(), rectangle).getInfo()[band]
        vmax = clipped_raster.reduceRegion(reducer = ee.Reducer.max(),
                                         geometry = rectangle,
                                         bestEfford=True).getInfo()[band]
    
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
    
    mapfig.fit_bounds(bounds=[(left_bottom_lat, left_bottom_lon), (right_top_lat, right_top_lon)])
    
    mapfig.add_child(rech_colormap)

# =============================================================================
# save map
# =============================================================================
    mapfig.save(outputfile)
