#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Image Kalman Filtering for MODIS-Landsat fusion
# Josh Gray, 2016
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Prelims
library(raster)
library(parallel)
library(rgdal)
library(dlm)
library(tools)
library(caTools)
library(RColorBrewer)
source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# main processing
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
# get the landsat data and dates
landsat_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
landsat_evi2_files <- dir(landsat_data_dir, pattern="evi2_landsat_overlap.tif", rec=T, full=T)
landsat_evi2_files <- landsat_evi2_files[order(GetLandsatDate(landsat_evi2_files))]
landsat_blue_files <- dir(landsat_data_dir, pattern="band1_landsat_overlap.tif", rec=T, full=T)
landsat_blue_files <- landsat_blue_files[order(GetLandsatDate(landsat_blue_files))]
landsat_green_files <- dir(landsat_data_dir, pattern="band2_landsat_overlap.tif", rec=T, full=T)
landsat_green_files <- landsat_green_files[order(GetLandsatDate(landsat_green_files))]
landsat_red_files <- dir(landsat_data_dir, pattern="band3_landsat_overlap.tif", rec=T, full=T)
landsat_red_files <- landsat_red_files[order(GetLandsatDate(landsat_red_files))]
landsat_nir_files <- dir(landsat_data_dir, pattern="band4_landsat_overlap.tif", rec=T, full=T)
landsat_nir_files <- landsat_nir_files[order(GetLandsatDate(landsat_nir_files))]
landsat_dates <- sort(GetLandsatDate(landsat_nir_files))

#-------------------------------------------------------------------------------
# get the modis data and dates
modis_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODISOVERLAP"
modis_evi2_files <- dir(modis_data_dir, pattern="evi2_modis_overlap.tif", full=T)
modis_dates <- as.Date(unlist(lapply(modis_evi2_files, GetModisDate)), origin="1970-1-1")
modis_evi2_files <- modis_evi2_files[order(modis_dates)]
modis_blue_files <- dir(modis_data_dir, pattern="blue_modis_overlap.tif", full=T)
modis_blue_files <- modis_blue_files[order(modis_dates)]
modis_green_files <- dir(modis_data_dir, pattern="green_modis_overlap.tif", full=T)
modis_green_files <- modis_green_files[order(modis_dates)]
modis_red_files <- dir(modis_data_dir, pattern="red_modis_overlap.tif", full=T)
modis_red_files <- modis_red_files[order(modis_dates)]
modis_nir_files <- dir(modis_data_dir, pattern="nir_modis_overlap.tif", full=T)
modis_nir_files <- modis_nir_files[order(modis_dates)]
modis_dates <- sort(modis_dates)
