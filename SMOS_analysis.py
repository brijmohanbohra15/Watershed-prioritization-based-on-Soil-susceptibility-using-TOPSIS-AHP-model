#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import sys


from netCDF4 import Dataset
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon



import os 
# choose directory where .nc file is located
os.chdir('E:\\Senitinel 1-2_SoilMoisture_FieldData-Updated\\_BrijMohanBohra\\SMOS\\filter_data\\filter_main_grid')

# folder_path where .ncfile is located
folder_path = os.getcwd()
#     extent a list = [(left_lat,lower_long),(left_lat,upper_long),(right_lat,upper_lat),(right_lat,lower_log)]
extent = [( 75.671970,30.358505), ( 75.671970,31.412555), ( 76.761606,31.412555), ( 76.761606, 30.358505)]




def smos_nc_shp_clip(folder_path,extent, extent_name):
#     import library glob,sys,os, netCDF4(Dataset),numpy as np, geopandas as gpd, matplotlib.pylot as plt, shapely.geometry(Point, Polygon)
#     crs = string in proj4 format for WGS-84 which is crs for smos file
#     extent a list = [(left_lat,lower_long),(left_lat,upper_long),(right_lat,upper_lat),(right_lat,lower_log)]
#     Ignore warning about missing/empty geometries
    import warnings
    warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)
    dir_list = os.listdir(folder_path)
    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' 
#   Clip geodataframe for your extent
    clip = gpd.GeoDataFrame()
    clip['geometry'] = None # create a new colum named geometry
    polygon = Polygon(extent) # Create a Shapely polygon from the coordinate-tuple list
    clip.loc[0, 'geometry'] = polygon # Insert the polygon into 'geometry' -column at index 0
    clip.loc[0, 'location'] = extent_name  # Add a new column and insert name of extent
    clip.crs = crs # Set the GeoDataFrame's coordinate system to WGS84 (i.e. epsg code 4326) or using gdf.crs
    
    for file in dir_list:
        
#         import netcdf file
        nc_file = Dataset(file)
        
#     extract variables
        lat = nc_file.variables['Latitude'][:]
        long = nc_file.variables['Longitude'][:]
        sm = nc_file.variables['Soil_Moisture'][:] 
        
#         make pandadataframe outof variables
        col = {'lats': lat,
        'longs': long,
        'soil_moisture': sm}
        nc_file_df = pd.DataFrame(col)
        nc_file_sub_df = nc_file_df.dropna()
        nc_file_sub_df.reset_index(drop= True, inplace=True)
        
#       make geodataframe out of dataframe using lat and lon column and crs file 
        gdf = gpd.GeoDataFrame(nc_file_sub_df, geometry=gpd.points_from_xy(nc_file_sub_df.longs, nc_file_sub_df.lats), crs=crs)
        print(gdf.shape)
        
# Clip the data using GeoPandas clip
        points_clip = gpd.clip(gdf, clip)
        if not points_clip.empty:
            # Determine the output path for the Shapefile
            outfp = file[:-2]+"shp"
            # Write the data into that Shapefile
            points_clip.to_file(outfp)         
        else:
            pass
           
 smos_nc_shp_clip(folder_path, extent, 'urna_extent')




