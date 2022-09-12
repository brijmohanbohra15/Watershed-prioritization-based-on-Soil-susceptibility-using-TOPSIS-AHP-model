
This repository covers the codes on my case study on Watershed prioritization based on Soil susceptibility using TOPSIS-AHP model.


## SMOS_anlaysis.py 
file contains function for clipping the SMOS netcdf file and interpolate that for given extent

## FSD.py
This file contains function for estimating fresh snow depth using Co-polar difference method.

## AHP.py
this file contains python script of Anlaytical Heirarchical Process ( MCDM technique ) to determine weight of  morphometric parameters; susceptible to soil erosion.
The output weightage of morphomatric parameter will be used in TOPSIS MCDM technique.

## TOPSIS.py
This file contains python script of TOPSIS (MCDM technique) to determine sub-watershed ranking based on susceptibility to soil erosion.

## sensitivityAnalysis_soilMoisture.py
This file contains fucntion to determine the coefficient of regression between soil moisture and backscatter coefficient with respect to different NDVI
ranges and available crop
