#!/usr/bin/env python
# coding: utf-8

# In[1]:
# =============================================================================
# Import libraries
# =============================================================================

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import gdal
import os
os.chdir(r"H:\BrijMohanBohra_wrd\UAV_SAR\PolSAR20170206")

# In[2]:
# =============================================================================
# Mean filetering
# =============================================================================
def meanFilter(img,kernel):
    output = img  
    for i in tqdm(range(kernel,img.shape[0]-kernel)):
        for j in range(kernel,img.shape[1]-kernel):
            block = img[i-kernel:i+kernel+1, j-kernel:j+kernel+1]
            #print(block)
            ans = []
            for z in block:
                ans.append(mean(z))
            m = mean(ans)
           # print(m)
            #print(i,j)
            # m = np.mean(block,dtype=np.float32)
            output[i][j] = m
    return output

def mean(img):
    l = len(img)
    lsum = 0
    for i in range(len(img)):
        if np.isnan(img[i]):
            lsum += 0
        else:
            lsum += img[i]
    meanl = lsum/l
    return meanl

# In[3]:
# =============================================================================
# Function for georefrencing 
# =============================================================================
    
def geoRefrence_img(fn, phi, array):
    # image name should end with .tif
    
    geot = phi.GetGeoTransform()
    proj = phi.GetProjection()
    driver_tiff = gdal.GetDriverByName("GTiff")
    drySnowDepth = driver_tiff.Create(fn, xsize=phi.RasterXSize, \
                                      ysize=phi.RasterYSize,bands=1,\
                                      eType = gdal.GDT_Float32)
    drySnowDepth.SetGeoTransform(geot)
    drySnowDepth.SetProjection(proj)
    drySnowDepth.GetRasterBand(1).WriteArray(array)
    
    
    # In[5]:
# =============================================================================
# Function for deopolarizaion factor
# =============================================================================
def get_depolarisation_factor(axial_ratio, shape):
    """
    Calculation three depolarisation factors
    :param axial_ratio: Axial ratio
    :param shape: Particle shape to consider
    :return: Tuple containing three depolarisation factors in x, y, and z directions
    """

    depolarisation_factorx = depolarisation_factory = depolarisation_factorz = 1 / 3.
    if axial_ratio == 1:
        return depolarisation_factorx, depolarisation_factory, depolarisation_factorz
    if shape == 'o':
        eccentricity = np.sqrt(axial_ratio ** 2 - 1)
        depolarisation_factorz = (1 + eccentricity ** 2) * (eccentricity - np.arctan(eccentricity)) / eccentricity ** 3
        depolarisation_factorx = depolarisation_factory = 0.5 * (1 - depolarisation_factorz)
    elif shape == 'p':
        eccentricity = np.sqrt(1 - axial_ratio ** 2)
        depolarisation_factorx = ((1 - eccentricity ** 2) *
                                  (np.log((1 + eccentricity) / (1 - eccentricity))
                                   - 2 * eccentricity)) / (2 * eccentricity ** 3)
        depolarisation_factory = depolarisation_factorz = 0.5 * (1 - depolarisation_factorx)
    return depolarisation_factorx, depolarisation_factory, depolarisation_factorz
# =============================================================================
# Function for effective permittivity
# =============================================================================
def get_effective_permittivity(fvol, depolarisation_factor):
    """
    Calculate effective permittivity
    :param fvol: Snow volume fraction
    :param depolarisation_factor: Depolarisation factor
    :return: Effective permittivity
    """

    eff_diff = EFF_ICE - EFF_AIR
    eff = EFF_AIR * (1 + fvol * eff_diff / (EFF_AIR + (1 - fvol) * depolarisation_factor * eff_diff))
    return eff

dep_fact = get_depolarisation_factor(1.5, 'o')

fvol = FRESH_SNOW_DENSITY/ ICE_DENSITY

eff_y = get_effective_permittivity(fvol, dep_fact[1])

eff_z = get_effective_permittivity(fvol, dep_fact[2])

eff_x = get_effective_permittivity(fvol, dep_fact[0])


# In[6]:

EFF_AIR = 1.00059
EFF_ICE = 3.179
ICE_DENSITY = 0.917  # gm/cc
FRESH_SNOW_DENSITY = 0.07  # gm/cc
WAVELENGTH = 24  # cm
NO_DATA_VALUE = -32768

# In[4]:
    
phi = gdal.Open('CPD_170206.tif')
phi_arr = phi.ReadAsArray()
filter_phase = meanFilter(phi_arr, 21)
filter_phase_float32 = filter_phase.astype('float32')
LIA = gdal.Open('LIA_170206.tif')
lia_arr = LIA.ReadAsArray()
lia_arr[lia_arr ==-10000] = np.nan

# In[7]:
n2h =eff_x
n2v = eff_y*(np.cos(lia_arr))**2 + eff_z*(np.sin(lia_arr))**2
n2h_new = np.zeros(n2v.shape)

for i in tqdm(range(n2h_new.shape[0])):
    for j in range(n2h_new.shape[1]):
        if np.isnan(n2v[i][j]):
            n2h_new[i][j]= np.nan
        else:
            n2h_new[i][j] = n2h

# In[8]:

mask = np.logical_and(filter_phase > 0, n2h_new**0.5>n2v**0.5)
mask_new = mask.astype('float32')


# In[9]:
tau = np.ones(n2v.shape)
tau = ((n2v- (np.sin(lia_arr))**2)**0.5- (n2h_new - (np.sin(lia_arr))**2)**0.5)*mask_new
tau[tau ==0] = np.nan

# In[10]:
drySnow = np.ones(n2v.shape)
drySnow = -1*WAVELENGTH*phi_arr/(4*np.pi*tau)
drySnow_filter = meanFilter(drySnow,65)

# In[3]:
# Create a raster
fn = 'phi_filter.tif'

# get the geotransform of the extisting raster
geot = phi.GetGeoTransform()
proj = phi.GetProjection()
get a gdal driver to use for raster creation
driver_tiff = gdal.GetDriverByName("GTiff")
# create a new raster dataset
print('starting raster creation')
new_ds = driver_tiff.Create(fn, xsize=10, ysize=10, bands=1, \
                           eType = gdal.GDT_Float32)
drySnowDepth = driver_tiff.Create(fn, xsize=phi.RasterXSize, ysize=phi.RasterYSize,bands=1, eType = gdal.GDT_Float32)
drySnowDepth2 = driver_tiff.Create(fn2, xsize=phi.RasterXSize, ysize=phi.RasterYSize,bands=1, eType = gdal.GDT_Float32)
SD = driver_tiff.Create(fn5, xsize=phi.RasterXSize, ysize=phi.RasterYSize,bands=1, eType = gdal.GDT_Float32)
                       
# set the geotransform
drySnowDepth2.SetGeoTransform(geot)
SD.SetGeoTransform(geot)
set the projection
drySnowDepth2.SetProjection(proj)
SD.SetProjection(proj)
print('rater created')

drySnowDepth.GetRasterBand(1).WriteArray(filter_phase_float32)
drySnowDepth2.GetRasterBand(1).WriteArray(filter_phase)
    
