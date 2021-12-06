#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import sys

# import fiona
# from fiona.crs import from_epsg
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon


# In[3]:


import os 
os.chdir('E:\\Senitinel 1-2_SoilMoisture_FieldData-Updated\\_BrijMohanBohra\\SMOS\\filter_data\\filter_main_grid')


# In[21]:


sm1 = Dataset('SM_OPER_MIR_SMUDP2_20210620T000734_20210620T010054_700_001_1.nc')
sm2 = Dataset('SM_OPER_MIR_SMUDP2_20210620T123814_20210620T133127_700_001_1.nc')
sm1


# In[5]:


lats1 = sm1.variables['Latitude'][:]
lats2 = sm2.variables['Latitude'][:]


# In[6]:


longs1 = sm1.variables['Longitude'][:]
longs2 = sm2.variables['Longitude'][:]


# In[7]:


datas1 = sm1.variables['Soil_Moisture'][:]
datas2 = sm2.variables['Soil_Moisture'][:]


# In[8]:


data1 = {'lats': lats1,
        'longs': longs1,
        'soil_moisture': datas1}
data2 = {'lats': lats2,
        'longs': longs2,
        'soil_moisture':datas2}


# In[12]:


df1 = pd.DataFrame(data2)
df2 = pd.DataFrame(data2)


# In[13]:


df1


# In[15]:


df2


# In[70]:


folder_path = os.getcwd()
folder_path
dir_list = os.listdir(folder_path)
dir_list
extent = [( 75.671970,30.358505), ( 75.671970,31.412555), ( 76.761606,31.412555), ( 76.761606, 30.358505)]
# for file in dir_list:
# #     nc_file = Dataset(file)
#     print(file[:-2]+"shp")
# clip = gpd.GeoDataFrame()
# clip['geometry'] = None # create a new colum named geometry
# polygon = Polygon(extent) # Create a Shapely polygon from the coordinate-tuple list
# clip.loc[0, 'geometry'] = polygon # Insert the polygon into 'geometry' -column at index 0
# clip.loc[0, 'location'] = 'extent_name'  # Add a new column and insert name of extent
# clip.crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
# print(clip.crs)
directory = 'clip_file'
path = os.path.join(folder_path, directory)
print(path)


# In[89]:



def smos_nc_shp_clip(folder_path,extent, extent_name):
#     import library glob,sys,os, netCDF4(Dataset),numpy as np, geopandas as gpd, matplotlib.pylot as plt, shapely.geometry(Point, Polygon)
# crs = string in proj4 format for WGS-84 which is crs for smos file
# extent a list = [(left_lat,lower_long),(left_lat,upper_long),(right_lat,upper_lat),(right_lat,lower_log)]
#  make new directory in folder name clip_file
#     directory = 'clip_file'
#     path = os.path.join(folder_path, directory)
#     os.mkdir(path) 
#     print(os.getcwd())
    # Ignore warning about missing/empty geometries
    import warnings
    warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)
    dir_list = os.listdir(folder_path)
    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' 
#     Clip geodataframe for your extent
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
#         make geodataframe out of dataframe using lat and lon column and crs file 
        gdf = gpd.GeoDataFrame(nc_file_sub_df, geometry=gpd.points_from_xy(nc_file_sub_df.longs, nc_file_sub_df.lats), crs=crs)
#         print(gdf.shape)
# Clip the data using GeoPandas clip
        points_clip = gpd.clip(gdf, clip)
        if not points_clip.empty:
            # Determine the output path for the Shapefile
            outfp = file[:-2]+"shp"
    # Write the data into that Shapefile
#             savefile = os.path.join(path,outfp)
            points_clip.to_file(outfp)         
        else:
            pass
    #do something

        
        
        
        


# In[90]:


smos_nc_shp_clip(folder_path, extent, 'urna_extent')


# In[ ]:





# In[ ]:




