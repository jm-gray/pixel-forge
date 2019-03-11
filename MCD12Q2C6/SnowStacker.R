#-------------------------------------------------------------------------------
C6SnowProcess <- function(x, tmp_download_dir, output_dir, myd11_qc_lut, overwrite=F, cleanup=T, nbar_scale_factor=1, lst_scale_factor=10, myd11a1=T, out_prefix="C6Snow"){
    # Final output 7-band stack: EVI2, NDSI, NDMI, A2 QA/QC, A2 Snow Flag, LST, LST Quality
    # return flags: 0=all processed, 1=SVI no LST, 2=LST no SVI, 3=nothing, 4=already exists and not overwritten

    # check for existence of output file
    out_file <- file.path(output_dir, paste(out_prefix, as.character(x[3]), as.character(x[1]), as.character(formatC(as.integer(x[2]), width=3, flag="0")), "tif", sep="."))
    if(file.exists(out_file) & !overwrite){
      return(4) # not processed because already exists and overwrite is F
    }
    
    # download the MCD43 and MYD11 files
    mcd43_files <- DownloadMCD43(x, out_root_dir=tmp_download_dir, overwrite=overwrite)
    if(myd11a1){
      myd11_file <- DownloadMYD11A1(x, out_root_dir=tmp_download_dir, overwrite=overwrite)
    }else{
      myd11_file <- DownloadMYD11A2(x, out_root_dir=tmp_download_dir, overwrite=overwrite)
    }

    # check if we got the MYD11 file
    if(length(myd11_file) == 0) myd11_file <- NA

    # case 1: we are missing one or both MCD43 and MYD11; don't process return "3" flag
    if((is.na(mcd43_files[1]) | is.na(mcd43_files[2])) & is.na(myd11_file)){
        return(3)
    }

    # case 2: we have both MCD43 but no MYD11; process SVI and QA, but no LST
    if(!is.na(mcd43_files[1]) & !is.na(mcd43_files[2]) & is.na(myd11_file)){
        evi_stack <- CalcEVI2(mcd43a4_file_path=mcd43_files[1], mcd43a2_file_path=mcd43_files[2], scale_factor=nbar_scale_factor, output_dir=NULL)
        ndsi <- CalcNDSI(mcd43a4_file_path=mcd43_files[1])
        ndmi <- CalcNDMI(mcd43a4_file_path=mcd43_files[1])
        blank_raster <- raster(evi_stack, 1)
        values(blank_raster) <- NA
        out_s <- stack(evi_stack[[1]], ndsi, ndmi, evi_stack[[2]], evi_stack[[3]], blank_raster, blank_raster)
        return_flag <- 1
    }

    # case 3: we have MYD11 but are missing one or both MCD43; calc LST but no SVI
    if((is.na(mcd43_files[1]) | is.na(mcd43_files[2])) & !is.na(myd11_file)){
        lst_stack <- ProcessMYD11QC(myd11_file, myd11_lut=myd11_qc_lut, lst_scale=lst_scale_factor, myd11a1=myd11a1)
        blank_raster <- raster(lst_stack, 1)
        values(blank_raster) <- NA
        out_s <- stack(blank_raster, blank_raster, blank_raster, blank_raster, blank_raster, lst_stack)
        return_flag <- 2
    }

    # case 4: we have both
    if(!is.na(mcd43_files[1]) & !is.na(mcd43_files[2]) & !is.na(myd11_file)){
        evi_stack <- CalcEVI2(mcd43a4_file_path=mcd43_files[1], mcd43a2_file_path=mcd43_files[2], scale_factor=nbar_scale_factor, output_dir=NULL)
        lst_stack <- ProcessMYD11QC(myd11_file, myd11_lut=myd11_qc_lut, lst_scale=lst_scale_factor, myd11a1=myd11a1)
        ndsi <- CalcNDSI(mcd43a4_file_path=mcd43_files[1])
        ndmi <- CalcNDMI(mcd43a4_file_path=mcd43_files[1])
        out_s <- stack(evi_stack[[1]], ndsi, ndmi, evi_stack[[2]], evi_stack[[3]], lst_stack)
        return_flag <- 0
    }

    # write the output to disk and return success flag
    writeRaster(out_s, file=out_file, datatype="INT2S")
    if(cleanup) system(paste("rm", paste(mcd43_files, collapse=" "), myd11_file))
    return(return_flag)
}

#-------------------------------------------------------------------------------
DownloadMCD43 <- function(x, out_root_dir, overwrite=F){
    # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
    out_dir <- file.path(out_root_dir, "MCD43A4", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
    existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    if(length(existing_mcd43a4_files) == 0 | overwrite){
        wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A4/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
        system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
        # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
        existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    }
    out_dir <- file.path(out_root_dir, "MCD43A2", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
    existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    if(length(existing_mcd43a2_files) == 0 | overwrite){
        wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A2/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
        system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
        # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
        existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    }
    return(c(existing_mcd43a4_files, existing_mcd43a2_files))
}

#-------------------------------------------------------------------------------
DownloadMYD11A2 <- function(x, out_root_dir, overwrite=F){
    # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
    out_dir <- file.path(out_root_dir, "MYD11A2", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
    existing_myd11a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    if(length(existing_myd11a2_files) == 0 | overwrite){
        wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MYD11A2/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
        system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
        # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
        existing_myd11a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    }
    return(existing_myd11a2_files)
}

#-------------------------------------------------------------------------------
DownloadMYD11A1 <- function(x, out_root_dir, overwrite=F){
    # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
    out_dir <- file.path(out_root_dir, "MYD11A1", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
    existing_myd11a1_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    if(length(existing_myd11a1_files) == 0 | overwrite){
        wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MYD11A1/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
        system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
        # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
        existing_myd11a1_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
    }
    return(existing_myd11a1_files)
}

#-------------------------------------------------------------------------------
DownloadMCD12Q1 <- function(x, out_root_dir, overwrite=F){
    # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
    # x form : c(year, tile)
    out_dir <- file.path(out_root_dir, "MCD12Q1", as.character(x[1]))
    existing_mcd12q1_files <- dir(out_dir, pattern=as.character(x[2]), full=T)
    if(length(existing_mcd12q1_files) == 0 | overwrite){
        wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[2]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD12Q1/', x[1], '/001', '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
        system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
        # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
        existing_mcd12q1_files <- dir(out_dir, pattern=as.character(x[2]), full=T)
    }
    return(existing_mcd12q1_files)
}

#-------------------------------------------------------------------------------
CalcNDSI <- function(mcd43a4_file_path, scale_factor=1, output_dir=NULL){
    # Computes NBAR-NDSI as: NDSI = B4 - B6 / B4 + B6
    # returns the raster or writes to disk depending on whether output_dir is NULL or not
    mcd43a4_file_path <- normalizePath(mcd43a4_file_path)
    b4_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=4, qual=F)
    b6_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=6, qual=F)

    b4_r <- raster(b4_sds)
    b6_r <- raster(b6_sds)

    ndsi <- ((b4_r * scale_factor) - (b6_r * scale_factor)) / ((b4_r * scale_factor) + (b6_r * scale_factor))
    ndsi <- round(ndsi * 1e4)
    plot(ndsi)

    if(!is.null(output_dir)){
        output_name <- file.path(output_dir, gsub(".*A([0-9]{7}).*(h[0-9]{2}v[0-9]{2}).*", "NDSI.\\1.\\2.tif", basename(mcd43a4_file_path)))
        writeRaster(ndsi, file=output_name, datatype="INT2S")
    }else{
        return(ndsi)
    }
}

#-------------------------------------------------------------------------------
CalcNDMI <- function(mcd43a4_file_path, scale_factor=1, output_dir=NULL){
    # Computes NBAR-NDSI as: NDSI = B4 - B6 / B4 + B6
    # returns the raster or writes to disk depending on whether output_dir is NULL or not
    mcd43a4_file_path <- normalizePath(mcd43a4_file_path)
    nir_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=2, qual=F)
    swir_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=6, qual=F)

    nir_r <- raster(nir_sds)
    swir_r <- raster(swir_sds)

    ndmi <- ((nir_r * scale_factor) - (swir_r * scale_factor)) / ((nir_r * scale_factor) + (swir_r * scale_factor))
    ndmi <- round(ndmi * 1e4)

    if(!is.null(output_dir)){
        output_name <- file.path(output_dir, gsub(".*A([0-9]{7}).*(h[0-9]{2}v[0-9]{2}).*", "NDMI.\\1.\\2.tif", basename(mcd43a4_file_path)))
        writeRaster(ndmi, file=output_name, datatype="INT2S")
    }else{
        return(ndmi)
    }
}

#-------------------------------------------------------------------------------
CalcEVI2 <- function(mcd43a4_file_path, mcd43a2_file_path, scale_factor=1, output_dir=NULL, out_prefix="C6EVI2"){
    # blue: 3, red: 1, nir: 2
    # Modify so returns EVI2, Qual Flag, Snow Flag
    blue_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=3, qual=F)
    blue_sds_qa <- GetMCD43A2SDSname(mcd43a2_file_path, band=3) # get the QA just from blue band
    red_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=1, qual=F)
    nir_sds <- GetMCD43A4SDSname(mcd43a4_file_path, band=2, qual=F)
    snow_sds <- GetMCD43A2SDSname(mcd43a2_file_path, band=3, snow=T) 

    # get RasterLayers
    blue_r <- raster(blue_sds)
    red_r <- raster(red_sds)
    nir_r <- raster(nir_sds)
    snow_r <- raster(snow_sds)
    blue_qa_r <- raster(blue_sds_qa)

    # calculate EVI2, EVI (if requested), and mask out snow values
    evi2 <- 2.5 * (((nir_r * scale_factor) - (red_r * scale_factor)) / ((nir_r * scale_factor) + (2.4 * (red_r * scale_factor)) + 1))
    # evi2[snow_r == 1] <- NA # this would fill snow values w/ NA, but we don't want to
    evi2 <- round(evi2 * 1e4)
    out_s <- stack(evi2, blue_qa_r, snow_r)

    # create the output stack and write to file or return
    if(!is.null(output_dir)){
        output_name <- file.path(output_dir, gsub(".*A([0-9]{7}).*(h[0-9]{2}v[0-9]{2}).*", paste(out_prefix, ".\\1.\\2.tif", sep=""), basename(mcd43a4_file_path)))
        writeRaster(out_s, file=output_name, datatype="INT2S")
    }else{
        return(out_s)
    }
}

#------------------------------------------------------------------------
GetMCD43A4SDSname <- function(mcd43a4_file_path, band, qual=F){
    if(qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_Band", band, sep=""), sep = ""))
    }else{
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
    }
}

#--------------------------------------------------------------------------------
GetMCD43A2SDSname <- function(mcd43a2_file_path, band, snow=F){
    # returns the SDS name for the Snow BRDF Albedo layer of an MCD43A2 file
    if(snow){
        # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:Snow_BRDF_Albedo
        return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
    }
    # HDF4_EOS:EOS_GRID:"MCD43A2.A2002002.h17v02.006.2016123222843.hdf":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band4
    return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band", band, sep=""), sep = ""))
}

#------------------------------------------------------------------------
GetMYD11SDSname <- function(myd11a2_file_path, day=T, qual=F, a1=F){
    if(a1){
        a1_a2_text="Daily"
    }else{
        a1_a2_text="8Day"
    }
    if(day & !qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", myd11a2_file_path, "\":MODIS_Grid_", a1_a2_text, "_1km_LST:LST_Day_1km", sep = ""))
    }else if(!day & !qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", myd11a2_file_path, "\":MODIS_Grid_", a1_a2_text, "_1km_LST:LST_Night_1km", sep = ""))
    }else if(day & qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", myd11a2_file_path, "\":MODIS_Grid_", a1_a2_text, "_1km_LST:QC_Day", sep = ""))
    }else if(!day & qual){
        return(paste("HDF4_EOS:EOS_GRID:\"", myd11a2_file_path, "\":MODIS_Grid_", a1_a2_text, "_1km_LST:QC_Night", sep = ""))
    }
}

#------------------------------------------------------------------------
ProcessMYD11QC <- function(myd11_file, myd11_lut, lst_scale=10, myd11a1=T){
    myd11_lst <- raster(GetMYD11SDSname(myd11_file, day=T, qual=F, a1=myd11a1))
    myd11_qc <- raster(GetMYD11SDSname(myd11_file, day=T, qual=T, a1=myd11a1))
    # process QC
    myd11_new_qc <- reclassify(myd11_qc, myd11_lut)
    # scale and restack
    myd11_s <- stack(round(myd11_lst * lst_scale), myd11_new_qc)
    # go to 500 m
    myd11_s <- disaggregate(myd11_s, fact=2)
    return(myd11_s)
}

#------------------------------------------------------------------------
bitFunc <- function(x, l=8) as.integer(unlist(lapply(as.character(intToBits(x)), function(y) strsplit(y, split="")[[1]][2])))[1:l]

#------------------------------------------------------------------------
MYD11QCtoInt <- function(x){
    # NA: LST not produced, 0: produced, LST error <= 1k, 1: produced, LST error <= 2k, 2: produced, LST error <= 3k, 3: produced, LST error > 3k
    bits <- bitFunc(x)
    mandatory_qa_flags <- paste(rev(bits[1:2]), collapse="")
    data_quality_flags <- paste(rev(bits[3:4]), collapse="")
    emis_error_flags <- paste(rev(bits[5:6]), collapse="")
    lst_error_flags <- paste(rev(bits[7:8]), collapse="")
    if(mandatory_qa_flags == "00"){
        return(0)
    }else if((mandatory_qa_flags == "10") | (mandatory_qa_flags == "11")){
        return(NA)
    }else if(lst_error_flags == "01"){
        return(1)
    }else if(lst_error_flags == "10"){
        return(2)
    }else if(lst_error_flags == "11"){
        return(3)
    }else{
        # captures cases where the only problem is in emissivity, but LST error is fine
        return(0)
    }
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(raster)
library(argparse)

arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-start", type="character") # start date as %Y-%m-%d
arg_parser$add_argument("-end", type="character") # end date as %Y-%m-%d
cmd_args <- arg_parser$parse_args()
# cmd_args <- arg_parser$parse_args(c("-tile","h12v03","-start","2013-1-1", "-end", "2013-12-31"))

# prelims
output_dir <- "/rsstu/users/j/jmgray2/SEAL/C6SNOW"
tmp_download_dir <- "/rsstu/users/j/jmgray2/SEAL"
log_file_dir <- "/rsstu/users/j/jmgray2/SEAL/C6SNOW/logs"

# create the MYD11 QC LUT
myd11_qc_lut <- as.matrix(data.frame(v=0:(2^8), qa=sapply(0:(2^8), MYD11QCtoInt)))

# create the full download data.frame
start_date <- as.Date(cmd_args$start)
end_date <- as.Date(cmd_args$end)
daily_dates <- seq.Date(start_date, end_date, by="day")
download_df <- data.frame(year=as.integer(strftime(daily_dates, format="%Y")), doy=as.integer(strftime(daily_dates, format="%j")), tile=cmd_args$tile)
download_df$tile <- as.character(download_df$tile)

process_flag <- apply(download_df, 1, C6SnowProcess, tmp_download_dir=tmp_download_dir, output_dir=output_dir, myd11_qc_lut=myd11_qc_lut)
# process_flag <- apply(download_df[1:3,], 1, C6SnowProcess, tmp_download_dir=tmp_download_dir, output_dir=output_dir, myd11_qc_lut=myd11_qc_lut)
log_file <- file.path(log_file_dir, paste("LogFile", cmd_args$tile, "Rdata", sep="."))
save(process_flag, download_df, file=log_file)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Rscript SnowStacker.R -tile h12v03 -start 2011-1-1 -end 2011-12-31 &
# Rscript SnowStacker.R -tile h12v03 -start 2012-1-1 -end 2012-12-31 &
# Rscript SnowStacker.R -tile h12v03 -start 2013-1-1 -end 2013-12-31 &
# Rscript SnowStacker.R -tile h12v03 -start 2014-1-1 -end 2014-12-31 &
