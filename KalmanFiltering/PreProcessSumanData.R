library(raster)

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
ScreenARDPIXELQA <- function(x, oli=FALSE){
    # parse PIXELQA value; returns 0 if: fill, cloud, cloud shadow, snow, water, med/high cloud conf,
    # or med/high confidence cirrus (oli=TRUE only). Otherwise, returns 1 
    # ScreenARDPIXELQA(c(1, 66, 68, 72, 80, 96, 130, 132, 136, 144, 160, 224)) # From Table 2-7 in ARD DFCB v6
    if(oli){
        as.integer(0 == bitwAnd(bitwShiftR(x, 0), 1) + bitwAnd(bitwShiftR(x, 2), 1) + bitwAnd(bitwShiftR(x, 3), 1) + bitwAnd(bitwShiftR(x, 4), 1) + bitwAnd(bitwShiftR(x, 5), 1) + bitwAnd(bitwShiftR(x, 7), 1)+ bitwAnd(bitwShiftR(x, 9), 1))
    }else{
        as.integer(0 == (bitwAnd(bitwShiftR(x, 0), 1) + bitwAnd(bitwShiftR(x, 2), 1) + bitwAnd(bitwShiftR(x, 3), 1) + bitwAnd(bitwShiftR(x, 4), 1) + bitwAnd(bitwShiftR(x, 5), 1) + bitwAnd(bitwShiftR(x, 7), 1)))
    }
}

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
ScreenARDRADSATQA <- function(x, oli=FALSE, bands=NULL){
    # returns 0 when the binary RADSATQA flag is on for ANY of the specified bands.
    # Defaults to all specral bands (1,2,3,4,5,7 for L4-7, 1:7 for L8)
    if(oli){
        if(is.null(bands)) bands <- 1:7
        as.integer(!any(bitwAnd(x, 2^bands)))
    }else{
        if(is.null(bands)) bands <- c(1:5, 7)
        as.integer(!any(bitwAnd(x, 2^bands)))
    }
}

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
GetARDMaskRaster <- function(pixelqa_file, radsatqa_file, bands=NULL, false_is_na=TRUE){
    # returns binary "use/dont't use" raster for masking ARD imagery
    # "1" only if not: fill, water, snow, cloud, cloud shadow, med/high cloud and cirrus conf
    # will be "0" otherwise, unless "false_is_na", in which case "don't use" will be NA
    # defaults to all bands (including 9-11 for Landsat)
    oli <- grepl("L(C|O)", basename(pixelqa_file)) # are these oli or oli+tirs ARD files?
    pixelqa_r <- raster(pixelqa_file)
    radsatqa_r <- raster(radsatqa_file)
    imagemask_r <- pixelqa_r
    pixelqa_mask_v <- ScreenARDPIXELQA(values(pixelqa_r), oli=oli)
    radsatqa_mask_v <- ScreenARDRADSATQA(values(radsatqa_r), oli=oli, bands=bands)
    values(imagemask_r) <- pixelqa_mask_v | radsatqa_mask_v
    if(false_is_na) imagemask_r[imagemask_r == 0] <- NA
    return(imagemask_r)
}

# #--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
# CropExample <- function(r, example_r, out_proj=NULL){
#     if(is.null(out_proj)) out_proj <- projection(example_r)
#     output_template_r <- projectExtent(example_r, out_proj)
#     crop_r <- crop(r, projectExtent(example_r, projection(r)), snap="out")
#     projectRaster()
#     return(crop_r)
#     # gdalwarp_params <- paste("-t_srs", out_proj, "-ts", ncol(example_r), nrow(example_r), "-te", paste(out_ex[c(1, 3, 2, 4)], collapse=" "), "-r cubic")
#     # gdal_cmd <- paste("gdalwarp", gdalwarp_params, in_file, out_file)
# }

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
CropResample <- function(r, example_r){
    # crops and resamples r to match example_r CRS, resolution, and extent
    crop_r <- crop(r, projectExtent(example_r, projection(r)))
    crop_proj_r <- projectRaster(crop_r, example_r)
    return(crop_proj_r)
}

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
PreprocessARD <- function(ard_dir, example_r, bands=NULL){
    # takes an ARD directory and constructs a multispectral image of "usable"
    # data that matches the geometry of example_r. If "bands" is NULL (default)
    # then a 6-band RasterStack: blue, green, red, nir, swir1, swir2 will be
    # returned. Otherwise, the returned RasterStack will be in order specified by "bands"

    # check if we're doing L8 or L4-7, and if bands are specified
    l8 <- grepl("L(C|O).*", basename(ard_dir))
    if(is.null(bands)){
        if(l8){
            # Landsat 8
            bands <- 2:7
        }else{
            # Landsat 4-7
            bands <- c(1:5, 7)
        }
    }

    # get the requested spectral bands and make a multispectral raster
    band_files <- dir(ard_dir, pattern=paste(".*SRB(", paste(bands, collapse="|"), ")", sep=""), full=T)
    band_nums <- as.integer(gsub(".*SRB([0-9]).*", "\\1", basename(band_files)))
    image_s <- stack(band_files[match(band_nums, bands)])
    
    # create the PIXELQA and RADSATQA mask and use it to screen the image
    pixelqa_file <- dir(ard_dir, pattern=".*PIXELQA.*.tif$", full=T)
    radsatqa_file <- dir(ard_dir, pattern=".*RADSATQA.*.tif$", full=T)
    imagemask_r <- GetARDMaskRaster(pixelqa_file = pixelqa_file, radsatqa_file = radsatqa_file, bands = bands)
    image_s <- image_s * imagemask_r

    # crop and reproject to match example_r
    image_crop_s <- CropResample(image_s, example_r)
    
    # eliminate values outside the valid range and apply scale factor
    image_crop_s[image_crop_s < 0 | image_crop_s > 10000] <- NA
    image_crop_s <- image_crop_s * 0.0001

    # return the final preprocessed image
    return(image_crop_s)
}

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
GetMCD43A4SDSname <- function(mcd43a4_file_path, band, qual=F){
  if(qual){
    return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_Band", band, sep=""), sep = ""))
  }else{
    return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
  }
}

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
GetMCD43A2SDSname <- function(mcd43a2_file_path, band, snow=F){
  # returns the SDS name for the Snow BRDF Albedo layer of an MCD43A2 file
  if(snow){
    # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:Snow_BRDF_Albedo
    return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
  }
  # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band4
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band", band, sep=""), sep = ""))
}


#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
# Prelims
# use the old extent to find new AEA (ARD native proj) extent and make example output raster
old_r <- raster(normalizePath("~/Google Drive File Stream/My Drive/Projects/KalmanFiltering/KF_fusion_data_new/landsat_new/LT50290312010329PAC01_sr_band2_landsat_overlap_sub.tif"))
example_ard_r <- raster(normalizePath("~/Desktop/ARD_Data/LC08_CU_016008_20180827_20190614_C01_V01_SR/LC08_CU_016008_20180827_20190614_C01_V01_SRB4.tif"))
out_proj <- projection(example_ard_r)
example_r <- crop(example_ard_r, projectExtent(old_r, crs=out_proj), snap="out")

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
# Landsat
ard_data_dir <- normalizePath("~/Desktop/ARD_Data")
ard_dirs <- dir(ard_data_dir, full=T)[!grepl(".tar$", dir(ard_data_dir))]

system.time(tmp <- PreprocessARD(ard_dirs[1], example_r = example_r))

# tmp_crop <- CropExample(tmp, r_landsat, projection(tmp))
system.time(tmp1_r <- PreprocessARD(ard_dirs[1], example_r = example_r))
system.time(tmp2_r <- PreprocessARD(ard_dirs[2], example_r = example_r))
system.time(tmp3_r <- PreprocessARD(ard_dirs[3], example_r = example_r))
layout(matrix(1:4, nrow=2))
par(mar=rep(0.1, 4), oma=rep(0, 4))
plotRGB(tmp1_r, 5, 4,3 , stretch="lin")
plotRGB(tmp2_r, 5, 4,3 , stretch="lin")
plotRGB(tmp3_r, 5, 4,3 , stretch="lin")

#--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+--*--+
# MODIS
a4_file <- normalizePath("~/Desktop/MCD43A4.A2018156.h10v04.006.2018166205545.hdf")
a2_file <- normalizePath("~/Desktop/MCD43A2.A2018156.h10v04.006.2018166205545.hdf")
modis_band_order <- c(3, 4, 1, 2, 6, 7)
modis_s <- stack(GetMCD43A4SDSname(a4_file, band=modis_band_order))

system.time(tmp <- CropResample(modis_s, example_r))
# plotRGB(modis_s, 2, 1, 3, stretch="lin")
plotRGB(tmp, 5, 4, 3, stretch="lin")

GetMCD43A2SDSname(a2_file, band=1:4)
GetMCD43A2SDSname(a2_file, band=1:4, snow=T)

r <- raster(GetMCD43A4SDSname(a4_file, band=1:4))


360/8

2 * 60 * tan(22.5 * (pi / 180))
0.70563
11/16
23/32
5/8
6/8