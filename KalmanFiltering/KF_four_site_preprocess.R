#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Download MODIS data and preprocess for 4 study site KF/STARFM comparison
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# The sites, Landsat path/rows and MODIS tiles:
# ---------------------------
# Central Valley, CA (Fresno/Biola)
# ---------------------------
# PR: 42, 35
# PR: 43, 34
# MODIS: h08v05
# ---------------------------
# Garden City, KS (Ulysses)
# ---------------------------
# PR: 30, 34
# PR: 31, 34
# MODIS: h09v05
# ---------------------------
# Adirondacks, NY
# ---------------------------
# PR: 14, 30
# PR: 15, 30
# MODIS: h12v04
# ---------------------------
# New York, NY
# ---------------------------
# PR: 13, 32
# PR: 14, 32
# MODIS: h12v04
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Prelims
library(parallel)
library(tools)
library(raster)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Functions
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
PlotOut <- function(i, files, dates, out_dir, prefix){
    s <- stack(files[i])
    out_jpeg <- file.path(out_dir, paste(prefix, formatC(i, width=3, flag="0"), ".jpg", sep=""))
    jpeg(file=out_jpeg, height=nrow(s), width=ncol(s), quality=100)
    plotRGB(s, 5, 4, 3, stretch="lin")
    text(x=mean(par()$usr[1:2]), y=par()$usr[4], label=dates[i], pos=1, cex=2, col="red")
    dev.off()    
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
bitFunc <- function(x, l=16) as.integer(unlist(lapply(as.character(intToBits(x)), function(y) strsplit(y, split="")[[1]][2])))[1:l]

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
remap_pixel_qa <- function(x, sensor="OLI"){
    # remaps pixel qa values
    # clear: 0
    # cloud: 1
    # shadow: 2
    # water: 3
    # snow: 4
    # fill: 5

    # handle the case of nodata
    if(is.na(x)) return(5) # map to fill

    bits <- bitFunc(x)
    if(sensor == "OLI"){
        # OLI
        if(bits[1] == 1){
            remap_qa <- 5 # fill
        }else if(bits[3] == 1){
            remap_qa <- 3 # water
        }else if(bits[4] == 1){
            remap_qa <- 2 # cloud shadow
        }else if(bits[5] == 1){
            remap_qa <- 4 # snow
        }else if(bits[6] == 1){
            remap_qa <- 1 # cloud
        }else if(any(paste(bits[10:9], collapse="") == c("11", "10"))){
            # map high/medium confidence cirrus to cloud
            remap_qa <- 1
        }else{
            remap_qa <- 0 # clear
        }
    }else{
        # TM, ETM+
        if(bits[1] == 1){
            remap_qa <- 5 # fill
        }else if(bits[3] == 1){
            remap_qa <- 3 # water
        }else if(bits[4] == 1){
            remap_qa <- 2 # cloud shadow
        }else if(bits[5] == 1){
            remap_qa <- 4 # snow
        }else if(bits[6] == 1){
            remap_qa <- 1 # cloud
        }else if(any(paste(bits[8:7], collapse="") == c("11", "10"))){
            # map medium and high cloud confidence to cloud
            remap_qa <- 1
        }else{
            remap_qa <- 0 # clear
        }
    }
    return(remap_qa)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
MakeTempDirectory <- function(output_root_dir){
  this_tmp_dir <- file.path(output_root_dir, paste(strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), sep="_"))
  while(dir.exists(this_tmp_dir)) this_tmp_dir <- file.path(output_root_dir, paste(strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), sep="_"))
  system(paste("mkdir", this_tmp_dir))
  return(this_tmp_dir)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
PreProcessARD <- function(tar_file, crop_ext, out_dir, out_suffix, gdal_program_path="/Library/Frameworks/GDAL.framework/Programs/"){
    # create the OLI LUT
    pixel_qa_lut_oli <- as.matrix(data.frame(v=0:(2^11), qa=sapply(0:(2^11), remap_pixel_qa)))
    pixel_qa_lut_oli <- rbind(pixel_qa_lut_oli, c(NA, 5)) # add a remap for NA values to be fill

    # create the ETM/TM LUT
    pixel_qa_lut_etm <- as.matrix(data.frame(v=0:(2^11), qa=sapply(0:(2^11), remap_pixel_qa, sensor="ETM")))
    pixel_qa_lut_etm <- rbind(pixel_qa_lut_oli, c(NA, 5)) # add a remap for NA values to be fill

    # untar
    this_tmp_dir <- MakeTempDirectory(out_dir)
    sensor <- gsub("(L[T|C|E][0-9]{2}).*", "\\1", basename(tar_file))
    if(sensor == "LC08"){
        untar_expression <- "'*SRB[2-7].tif' '*PIXELQA.tif' '*.xml'"
        pixel_qa_lut <- pixel_qa_lut_oli
    }else{
        untar_expression <- "'*SRB[1-5,7].tif' '*PIXELQA.tif' '*.xml'"
        pixel_qa_lut <- pixel_qa_lut_etm
    }
    tar_cmd <- paste("tar -xvf", tar_file, "-C", this_tmp_dir, untar_expression)
    system(tar_cmd)

    # gather all files
    multispec_files <- dir(this_tmp_dir, pattern="SRB.*tif$", full=T)
    pixelqa_file <- dir(this_tmp_dir, pattern="PIXELQA.tif$", full=T)
    xml_file <- dir(this_tmp_dir, pattern="xml$", full=T)
    
    # crop and remap PIXELQA
    pixelqa_r <- raster(pixelqa_file)
    pixelqa_r <- crop(pixelqa_r, crop_ext)
    pixelqa_r <- reclassify(pixelqa_r, pixel_qa_lut)
    # check if all fill
    if(minValue(pixelqa_r) != 0){
        # DO NOTHING
        print("All fill, not writing output")
        return(NA)
    }
    out_pixelqa_file <- ""
    out_pixelqa_file <- file.path(out_dir, paste(file_path_sans_ext(basename(pixelqa_file)), out_suffix, ".envi", sep=""))
    writeRaster(pixelqa_r, file=out_pixelqa_file, datatype="INT2S", NAflag=32767, format="ENVI")

    # crop multispec files, screen with pixelqa and write to disk
    out_files <- c()
    for(this_file in multispec_files){
        # this_file <- multispec_files[1]
        out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(this_file)), out_suffix, ".envi", sep=""))
        r <- raster(this_file)
        r_crop <- crop(r, crop_ext)
        r_crop[pixelqa_r != 0] <- NA
        writeRaster(r_crop, file=out_file, datatype="INT2S", NAflag=32767, format="ENVI")
        out_files <- c(out_files, out_file)
    }
    out_files <- c(out_files, out_pixelqa_file)

    # create VRT for multispectral bands
    out_vrt_file <- file.path(out_dir, paste(gsub("(.*)_[A-Z0-9]{4}.*tif", "\\1", basename(this_file)), out_suffix, ".vrt", sep=""))
    gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), "-separate", out_vrt_file, paste(out_files, collapse=" "))
    system(gdalvrt_cmd)

    # remove temporary directory
    system(paste("rm -r", this_tmp_dir))

    return(out_vrt_file)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
DownloadMCD43 <- function(x, out_root_dir, overwrite=F, a4=T){
    # function to download MCD43A4 (a4=T) and MCD43A2 (a4=F) data
    if(a4){
        # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
        out_dir <- file.path(out_root_dir, "MCD43A4", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
        existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
        if(length(existing_mcd43a4_files) == 0 | overwrite){
            wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A4/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
            # system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
            system(wget_cmd)
            # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
            existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
            return(existing_mcd43a4_files)
        }
    }else{
        out_dir <- file.path(out_root_dir, "MCD43A2", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
        existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
        if(length(existing_mcd43a2_files) == 0 | overwrite){
            wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A2/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
            # system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
            system(wget_cmd)
            # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
            existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
        }
        return(existing_mcd43a2_files)
    }
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
CheckMCD43 <- function(x, out_root_dir, a4=T){
    # check if data are successfully downloaded
    if(a4){
        out_dir <- file.path(out_root_dir, "MCD43A4", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
        existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
        if(length(existing_mcd43a4_files) == 0){
            return(FALSE)
        }else{
            return(TRUE)
        }
    }else{
        out_dir <- file.path(out_root_dir, "MCD43A2", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
        existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
        if(length(existing_mcd43a2_files) == 0){
            return(FALSE)
        }else{
            return(TRUE)
        }
    }
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
PreProcessModis <- function(the_date, the_tile, a4dir, a2dir, out_dir, ex_r, out_suffix, gdal_program_path="/Library/Frameworks/GDAL.framework/Programs/"){
    this_tmp_dir <- MakeTempDirectory(out_dir)
    a4_file <- dir(file.path(a4dir, strftime(the_date, format="%Y"), strftime(the_date, format="%j")), pattern=paste(".*", tile, ".*", "hdf$", sep=""), full=T)
    a2_file <- dir(file.path(a2dir, strftime(the_date, format="%Y"), strftime(the_date, format="%j")), pattern=paste(".*", tile, ".*", "hdf$", sep=""), full=T)

    # reproject, resample, and crop the A2 snow and QA layers
    tmp_snow_file <- file.path(this_tmp_dir, paste(file_path_sans_ext(basename(a2_file)), "_snow_", ".tiff", sep=""))
    gdalwarp_cmd_snow <- paste(gdal_program_path, "gdalwarp -t_srs '", projection(ex_r),"' -te ", paste(extent(ex_r)[c(1, 3, 2, 4)], collapse=" "), " -ts ", ncol(ex_r), " ", nrow(ex_r), " ", GetMCD43A2SDSname(a2_file, band=3, snow=T), " ", tmp_snow_file, sep="")
    system(gdalwarp_cmd_snow)
    tmp_qa_file <- file.path(this_tmp_dir, paste(file_path_sans_ext(basename(a2_file)), "_qa_", ".tiff", sep=""))
    gdalwarp_cmd_qa <- paste(gdal_program_path, "gdalwarp -t_srs '", projection(ex_r),"' -te ", paste(extent(ex_r)[c(1, 3, 2, 4)], collapse=" "), " -ts ", ncol(ex_r), " ", nrow(ex_r), " ", GetMCD43A2SDSname(a2_file, band=3), " ", tmp_qa_file, sep="")
    system(gdalwarp_cmd_qa)
    r_snow <- raster(tmp_snow_file)
    r_qa <- raster(tmp_qa_file)
    r_mask <- r_snow
    values(r_mask) <- 0
    r_mask[r_snow != 0 | r_qa == 4] <- 1
    qa_out_file <- file.path(out_dir, paste("MCD43A4_", the_tile, "_", strftime(the_date, "%Y"), "_", strftime(the_date, "%j"), "_band_qa_", out_suffix, sep=""))
    writeRaster(r_mask, file=qa_out_file, datatype="INT2S", NAflag=32767, format="ENVI")

    # reproject, resample, and crop each A4 multispectral file, apply QA/snow mask, and write to disk
    band_names <- c("red", "nir", "blue", "green", "5", "swir1", "swir2")
    for(band in c(3, 4, 1, 2, 6, 7)){
        tmp_out_file <- file.path(this_tmp_dir, paste(file_path_sans_ext(basename(a4_file)), "_band_", band, "_", band_names[band], "_", out_suffix, ".tif", sep=""))
        gdalwarp_cmd <- paste(gdal_program_path, "gdalwarp -t_srs '", projection(ex_r),"' -te ", paste(extent(ex_r)[c(1, 3, 2, 4)], collapse=" "), " -ts ", ncol(ex_r), " ", nrow(ex_r), " ", GetMCD43A4SDSname(a4_file, band=band), " ", tmp_out_file, sep="")
        system(gdalwarp_cmd)
        r <- raster(tmp_out_file)
        r[r_mask != 0] <- NA
        out_file <- file.path(out_dir, paste("MCD43A4_", the_tile, "_", strftime(the_date, "%Y"), "_", strftime(the_date, "%j"), "_band_", band, "_", band_names[band], "_", out_suffix, sep=""))
        writeRaster(r, file=out_file, datatype="INT2S", NAflag=32767, format="ENVI")
        i <- i + 1
    }
    
    # create a VRT
    # out_files <- dir(out_dir, pattern=paste(".*", the_tile, "_", strftime(the_date, "%Y"), "_", strftime(the_date, "%j"), ".*envi$", sep=""), full=T)
    # out_vrt_file <- file.path(out_dir, paste("MCD43A4_", the_tile, "_", strftime(the_date, "%Y"), "_", strftime(the_date, "%j"), "_", out_suffix, ".vrt", sep=""))
    # gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), "-separate", out_vrt_file, paste(out_files, collapse=" "))
    # system(gdalvrt_cmd)
    
    # cleanup
    system(paste("rm -r", this_tmp_dir))
    return(out_vrt_file)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
WriteVRT <- function(the_date, the_tile, out_dir, out_suffix, gdal_program_path="/Library/Frameworks/GDAL.framework/Programs/"){
    out_files <- c()
    for(the_band in c(3, 4, 1, 2, 6, 7)){
        out_files <- c(out_files, (file.path(out_dir, paste(paste("MCD43A4", the_tile, strftime(the_date, format="%Y"), strftime(the_date, format="%j"), "band", the_band, band_names[the_band], out_suffix, sep="_"), ".envi", sep=""))))
    }
    out_files <- c(out_files, file.path(out_dir, paste(paste("MCD43A4", the_tile, strftime(the_date, format="%Y"), strftime(the_date, format="%j"), "band_qa", out_suffix, sep="_"), ".envi", sep="")))
    out_vrt_file <- file.path(out_dir, paste("MCD43A4_", the_tile, "_", strftime(the_date, "%Y"), "_", strftime(the_date, "%j"), "_", out_suffix, ".vrt", sep=""))
    gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), "-separate", out_vrt_file, paste(out_files, collapse=" "))
    system(gdalvrt_cmd)
}


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
GetMCD43A4SDSname <- function(mcd43a4_file_path, band, qual=F){
    if(qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_Band", band, sep=""), sep = ""))
    }else{
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
    }
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
GetMCD43A2SDSname <- function(mcd43a2_file_path, band, snow=F){
    # returns the SDS name for the Snow BRDF Albedo layer of an MCD43A2 file
    if(snow){
        # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:Snow_BRDF_Albedo
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
    }
    # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band4
    return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band", band, sep=""), sep = ""))
}


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Do the downloading
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
tiles <- c("h08v05", "h09v05", "h12v04")
start_date <- as.Date("2018-1-1")
end_date <- as.Date("2018-12-31")
daily_dates <- seq.Date(start_date, end_date, by="day")
download_df <- expand.grid(date=daily_dates, tile=tiles)
download_df$tile <- as.character(download_df$tile)
download_df$doy <- as.integer(strftime(download_df$date, format="%j"))
download_df$year <- as.integer(strftime(download_df$date, format="%Y"))
download_df <- download_df[,c(4,3,2)]

output_dir <- "/Volumes/users/j/jmgray2/SEAL/MCD43C6New"

trash <- apply(download_df, 1, DownloadMCD43, out_root_dir=output_dir, a4=T)
trash <- apply(download_df, 1, DownloadMCD43, out_root_dir=output_dir, a4=F)

a4_success <- apply(download_df, 1, CheckMCD43, out_root_dir=output_dir, a4=T)
a2_success <- apply(download_df, 1, CheckMCD43, out_root_dir=output_dir, a4=F)

trash <- apply(download_df[!a4_success,], 1, DownloadMCD43, out_root_dir=output_dir, a4=T)
trash <- apply(download_df[!a2_success,], 1, DownloadMCD43, out_root_dir=output_dir, a4=F)


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Define the subsets and cut the ARD data
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# 1) Binary output, 2-byte integer, scaled by 10000
# Landsat  MODIS
#    1       3
#    2       4
#    3       1
#    4       2
#    5       6
#    7       7

library(raster)
library(rgdal)

landsat_pathrow_shp_file <- normalizePath("~/Google Drive/DataSets/wrs2_descending/wrs2_descending.shp")
landsat_pathrow_shp <- shapefile(landsat_pathrow_shp_file)
modis_grid_shp_file <- normalizePath("/Users/jmgray2/Google Drive/DataSets/modis_grid/modis_grid_wgs84/modis_sinusoidal_grid_world.shp")
modis_grid_shp <- shapefile(modis_grid_shp_file)

# central valley
cv_landsat <- landsat_pathrow_shp[(landsat_pathrow_shp$PATH == 42 & landsat_pathrow_shp$ROW == 35) | (landsat_pathrow_shp$PATH == 43 & landsat_pathrow_shp$ROW == 34), ]
cv_modis <- modis_grid_shp[modis_grid_shp$h == 8 & modis_grid_shp$v == 5,]
writeOGR(cv_landsat, "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp", "cv_landsat", driver="ESRI Shapefile")
writeOGR(cv_modis, "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp", "cv_modis", driver="ESRI Shapefile")
# ogr2ogr -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' cv_landsat_aea cv_landsat.shp
# ogr2ogr -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' cv_modis_aea cv_modis.shp
cv_landsat_aea <- shapefile("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp/cv_landsat_aea/cv_landsat.shp")
cv_modis_aea <- shapefile("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp/cv_modis_aea/cv_modis.shp")
cv_ard_ex_s <- stack(dir("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/example_ARD/central_valley/LC08_CU_003010_20181126_20181216_C01_V01_SR/", pattern=".*SRB.*tif", full=T))
cut_size <- 22e3
cv_cut_ex <- extent(-2111300, -2111300 + cut_size, 1800680 - cut_size, 1800680)
plot(cv_modis_aea, lty=2)
plot(cv_landsat_aea, add=T, lty=1)
plot(cv_landsat_aea, lty=1)
plotRGB(cv_ard_ex_s, 5, 4, 3, stretch="lin", add=T)
plot(cv_cut_ex, add=T, col="yellow")
plot(cv_landsat_aea, lty=1, add=T, border="green")
out_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard/"
tar_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/"
tar_files <- dir(tar_dir, pattern="*SR.tar$", full=T, rec=T)
i <- 1
for(tar_file in tar_files){
    print(paste(i, "of", length(tar_files), ":", basename(tar_file)))
    the_vrt_file <- PreProcessARD(normalizePath(tar_file), crop_ext=cv_cut_ex, out_dir=out_dir, out_suffix="_cv_sub")
    i <- i + 1
}
# make a movie
cv_ard_files <- dir(out_dir, pattern=".vrt$", full=T)
cv_ard_sensors <- gsub("(L[T|C|E][0-9]{2}).*", "\\1", basename(cv_ard_files))
cv_ard_dates <- as.Date(gsub(".*_[0-9]{6}_([0-9]{8})_[0-9]{8}.*", "\\1", basename(cv_ard_files)), format="%Y%m%d")
cv_ard_files <- cv_ard_files[order(cv_ard_dates)]
cv_ard_sensors <- cv_ard_sensors[order(cv_ard_dates)]
cv_ard_dates <- cv_ard_dates[order(cv_ard_dates)]
trash <- lapply(1:length(cv_ard_files), PlotOut, files=cv_ard_files, dates=cv_ard_dates, out_dir=jpeg_dir, prefix="CV_ard_")
# ffmpeg -r 10 -i CV_ard_%03d.jpg -vcodec libx264 -y -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -an cv_ard.mp4 


# garden city
gc_landsat <- landsat_pathrow_shp[(landsat_pathrow_shp$PATH == 30 | landsat_pathrow_shp$PATH == 31) & landsat_pathrow_shp$ROW == 34, ]
gc_modis <- modis_grid_shp[(modis_grid_shp$h == 9 | modis_grid_shp$h == 10) & modis_grid_shp$v == 5,]
writeOGR(gc_landsat, "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp", "gc_landsat", driver="ESRI Shapefile")
writeOGR(gc_modis, "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp", "gc_modis", driver="ESRI Shapefile")
# ogr2ogr -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' gc_landsat_aea gc_landsat.shp
# ogr2ogr -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' gc_modis_aea gc_modis.shp
gc_landsat_aea <- shapefile("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp/gc_landsat_aea/gc_landsat.shp")
gc_modis_aea <- shapefile("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/study_shp/gc_modis_aea/gc_modis.shp")
gc_ard_ex_s <- stack(dir("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/example_ARD/garden_city/LC08_CU_013011_20181208_20190106_C01_V01_SR", , pattern=".*SRB.*tif", full=T))
plot(gc_modis_aea, lty=2)
plot(gc_landsat_aea, add=T, lty=1)
plot(gc_landsat_aea, lty=1)
plotRGB(gc_ard_ex_s, 5, 4, 3, stretch="lin", add=T)

# adirondacks
ad_landsat <- landsat_pathrow_shp[(landsat_pathrow_shp$PATH == 14 | landsat_pathrow_shp$PATH == 15) & landsat_pathrow_shp$ROW == 30, ]
ad_modis <- modis_grid_shp[modis_grid_shp$h == 12 & modis_grid_shp$v == 4,]
plot(ad_modis, lty=2)
plot(ad_landsat, add=T, lty=1)

# nyc
ny_landsat <- landsat_pathrow_shp[(landsat_pathrow_shp$PATH == 13 | landsat_pathrow_shp$PATH == 14) & landsat_pathrow_shp$ROW == 32, ]
ny_modis <- modis_grid_shp[modis_grid_shp$h == 12 & (modis_grid_shp$v == 4 | modis_grid_shp$v == 5),]
plot(ny_modis, lty=2)
plot(ny_landsat, add=T, lty=1)

##########################
site <- "central_valley"
root_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat"
in_files <- dir(file.path(root_dir, site), pattern="*SR.tar$", full=T, rec=T)
sensor <- gsub("(L[A-Z][0-9]{2})_.*", "\\1", basename(in_files))
ard_tile <- gsub("L.*_([0-9]{6})_.*", "\\1", basename(in_files))
dates <- as.Date(gsub("L.*_([0-9]{8})_[0-9]{6}.*", "\\1", basename(in_files)), format="%Y%m%d")

library(viridis)
plot_cols <- magma(5)[c(2,4)]
plot(dates, rep(1, length(dates)), pch=16, col=plot_cols[as.integer(as.factor(sensor))], xlim=range(dates))


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Preprocess MODIS
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
a4dir <- "/Volumes/users/j/jmgray2/SEAL/MCD43A4"
a2dir <- "/Volumes/users/j/jmgray2/SEAL/MCD43A2"

#--------------------------------
# Central Valley
out_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43"
ex_r <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard/LC08_CU_003010_20180315_20181129_C01_V01_cv_sub.vrt")
the_tile <- "h08v05"
out_suffix <- "cv_sub"
all_dates <- seq.Date(as.Date("2018-1-1"), as.Date("2018-12-31"), by="day")

i <- 1
for(the_date in all_dates){
    print(paste("Doing:", the_date, "|", i, "of", length(all_dates)))
    print("-----------------------------------------")
    trash_vrt <- PreProcessModis(the_date=as.Date(the_date, origin="1970-1-1"), the_tile=the_tile, a4dir=a4dir, a2dir=a2dir, out_dir=out_dir, ex_r=ex_r, out_suffix=out_suffix)
    i <- i + 1
}

# write the VRTs (because we flubbed it the first time...)
i <- 1
for(the_date in all_dates){
    print(paste("Doing:", the_date, "|", i, "of", length(all_dates)))
    print("-----------------------------------------")
    trash_vrt <- WriteVRT(the_date=as.Date(the_date, origin="1970-1-1"), the_tile=the_tile, out_dir=out_dir, out_suffix=out_suffix)
    i <- i + 1
}

#--------------------------------
# Garden City
out_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/garden_city/gc_mcd43"
ex_r <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard/LC08_CU_003010_20180315_20181129_C01_V01_cv_sub.vrt")
the_tile <- "h08v05"
out_suffix <- "cv_sub"
all_dates <- seq.Date(as.Date("2018-1-1"), as.Date("2018-12-31"), by="day")

i <- 1
for(the_date in all_dates){
    print(paste("Doing:", the_date, "|", i, "of", length(all_dates)))
    print("-----------------------------------------")
    trash_vrt <- PreProcessModis(the_date=as.Date(the_date, origin="1970-1-1"), the_tile=the_tile, a4dir=a4dir, a2dir=a2dir, out_dir=out_dir, ex_r=ex_r, out_suffix=out_suffix)
    i <- i + 1
}

# write the VRTs (because we flubbed it the first time...)
i <- 1
for(the_date in all_dates){
    print(paste("Doing:", the_date, "|", i, "of", length(all_dates)))
    print("-----------------------------------------")
    trash_vrt <- WriteVRT(the_date=as.Date(the_date, origin="1970-1-1"), the_tile=the_tile, out_dir=out_dir, out_suffix=out_suffix)
    i <- i + 1
}
