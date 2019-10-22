ard_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard"
mcd43_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43"

ard_vrts <- dir(ard_dir, pattern=".*vrt$", full=T)
ard_dates <- as.Date(gsub(".*_[0-9]{6}_([0-9]{8})_[0-9]{8}_.*", "\\1", basename(ard_vrts)), format="%Y%m%d")
s <- stack(ard_vrts[5])
plotRGB(s, 4, 3, 2, stretch="lin")

ard_dates[3]
ard_dates[5]
strftime(ard_dates[5], format="%Y%j")
ard_vrts[5]

/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard/LC08_CU_003010_20180126_20181203_C01_V01_cv_sub.vrt

Input ARD:
LC08_CU_003010_20180126_20181203_C01_V01_SRB2_cv_sub.envi
LC08_CU_003010_20180211_20181130_C01_V01_SRB4_cv_sub.envi
Match MCD:
MCD43A4_h08v05_2018_026_band_2_nir_cv_sub.envi
MCD43A4_h08v05_2018_042_band_2_nir_cv_sub.envi

Pred ARD: LC08_CU_003010_20180202_20181203_C01_V01_cv_sub
Pred MCD: MCD43A4_h08v05_2018_033_band_2_nir_cv_sub.envi

Input date: 1/26/2018
Pred date: 2/2/2018
strftime(as.Date("2018-1-26"), format="%Y-%j")
strftime(as.Date("2018-2-2"), format="%Y-%j")

# r_starfm <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/starfm/starfm_2018_033_nir.bin")
r_starfm <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/starfm/starfm_2018_033_nir_2match.bin")
r_actual <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard/LC08_CU_003010_20180202_20181203_C01_V01_SRB4_cv_sub.envi")
r_mod1 <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43/MCD43A4_h08v05_2018_026_band_2_nir_cv_sub.envi")
r_mod2 <- raster("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43/MCD43A4_h08v05_2018_033_band_2_nir_cv_sub.envi")
extent(r_starfm) <- extent(r_actual)
projection(r_starfm) <- projection(r_actual)
NAvalue(r_starfm) <- 32767
NAvalue(r_actual) <- 32767
NAvalue(r_mod1) <- 32767
NAvalue(r_mod2) <- 32767
r_diff <- r_actual - r_starfm
qs_starfm <- quantile(r_starfm, c(0, 0.02, 0.98, 1), na.rm=T)
qs_actual <- quantile(r_actual, c(0, 0.02, 0.98, 1), na.rm=T)
r_landsat_breaks <- c(min(qs_starfm[1], qs_actual[1]), seq(min(qs_starfm[2], qs_actual[2]), max(qs_starfm[3], qs_actual[3]), len=254), max(qs_starfm[4], qs_actual[4]))
qs_diff <- quantile(r_diff, c(0, 0.02, 0.98, 1), na.rm=T)
diff_breaks <- c(qs_diff[1], seq(-1 * max(abs(qs_diff[2:3])), max(abs(qs_diff[2:3])), len=254), qs_diff[4])
diff_pal <- colorRampPalette(brewer.pal(9, "RdBu"))
qs_mod1 <- quantile(r_mod1, c(0, 0.02, 0.98, 1))
qs_mod2 <- quantile(r_mod2, c(0, 0.02, 0.98, 1))
r_modis_breaks <- c(min(qs_mod1[1], qs_mod2[1]), seq(min(qs_mod1[2], qs_mod2[2]), max(qs_mod1[3], qs_mod2[3]), len=254), max(qs_mod1[4], qs_mod2[4]))

layout(matrix(1:6, nrow=2, byrow=T))
plot(r_mod1, col=plasma(255), maxpixels=ncell(r_mod1), breaks=r_modis_breaks)
title("MCD43 Match Day")
plot(r_mod2, col=plasma(255), maxpixels=ncell(r_mod1), breaks=r_modis_breaks)
title("MCD43 Pred Day")
hist(values(r_diff), xlab="ARD - StarFM", ylab="Count", col="grey", border="white", main="ARD - STARFM")
box()
plot(r_actual, col=plasma(255), maxpixels=ncell(r_actual), breaks=r_landsat_breaks)
title("ARD Pred Day")
plot(r_starfm, col=plasma(255), maxpixels=ncell(r_starfm), breaks=r_landsat_breaks)
title("StarFM Pred Day")
plot(r_diff, breaks=diff_breaks, col=diff_pal(255), maxpixels=ncell(r_diff))
title("ARD - STARFM")


landsat_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard"
modis_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43"
LC08_CU_003010_20180211_20181130_C01_V01_SRB4_cv_sub.envi
MCD43A4_h08v05_2018_026_band_2_nir_cv_sub.envi

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
MakeStarFM <- function(match_date_1, match_date_2=NA, pred_date, landsat_dir, modis_dir, out_dir, out_prefix="starfm", run_cmd=FALSE, use_spatial=TRUE, landsat_uncertainty=50, modis_uncertainty=50, max_search_dist=NA, num_slice=40, starfm_path="/Users/jmgray2/Google\ Drive/Projects/KalmanFiltering/STARFM/StarFM/source/StarFM.exe"){
    # creates a StarFM fun file and optionally executes it (run_cmd=TRUE)
    # for all 6 Landsat multispectral bands (blue, green, red, nir, swir1, swir2)
    # assumes that NA flags are 32767 for both Landsat and MODIS
    # assumes that scale factor is 10000
    # supports 1 or 2 input pairs, and spatial, uncertainty, and pure slice options
    band_names <- c("blue", "green", "red", "nir", "swir1", "swir2")
    
    # get MODIS and Landsat input files for first match date
    landsat_match_date_1_files <- dir(landsat_dir, pattern=paste(".*", strftime(match_date_1, format="%Y%m%d"), "_[0-9]{8}_.*SRB[0-9]{1}.*envi$", sep=""), full=T)
    if(length(landsat_match_date_1_files) != length(band_names)){
        print(paste("Did not find", length(band_names), "Landsat files for", match_date_1, "in", landsat_dir, "|| Aborting"))
        return(NA)
    }
    landsat_match_date_1_bands <- as.integer(gsub(".*_SRB([0-9]{1})_.*", "\\1", basename(landsat_match_date_1_files)))
    landsat_match_date_1_sensor <- gsub("L[C|E|T]([0-9]{2}).*", "\\1", basename(landsat_match_date_1_files[1]))
    if(landsat_match_date_1_sensor == "08"){
        landsat_match_date_1_band_list <- list("blue"=2, "green"=3, "red"=4, "nir"=5, "swir1"=6, "swir2"=7)
    }else{
        landsat_match_date_1_band_list <- list("blue"=1, "green"=2, "red"=3, "nir"=4, "swir1"=5, "swir2"=7)
    }
    modis_band_list <- list("blue"=3, "green"=4, "red"=1, "nir"=2, "swir1"=6, "swir2"=7)
    modis_match_date_1_files <- dir(modis_dir, pattern=paste(".*", strftime(match_date_1, format="%Y_%j"), "_band_[0-9]{1}.*envi$", sep=""), full=T)
    if(length(modis_match_date_1_files) != length(band_names)){
        print(paste("Did not find", length(band_names), "MODIS files for", match_date_1, "in", modis_dir, "|| Aborting"))
        return(NA)
    }
    modis_match_date_1_bands <- as.integer(gsub(".*_band_([0-9]{1})_.*", "\\1", basename(modis_match_date_1_files)))
    
    # get MODIS and Landsat input files for second match date, if provided
    if(!is.na(match_date_2)){
        landsat_match_date_2_files <- dir(landsat_dir, pattern=paste(".*", strftime(match_date_2, format="%Y%m%d"), "_[0-9]{8}_.*SRB[0-9]{1}.*envi$", sep=""), full=T)
        if(length(landsat_match_date_2_files) != length(band_names)){
            print(paste("Did not find", length(band_names), "Landsat files for", match_date_2, "in", landsat_dir, "|| Aborting"))
            return(NA)
        }
        landsat_match_date_2_bands <- as.integer(gsub(".*_SRB([0-9]{1})_.*", "\\1", basename(landsat_match_date_2_files)))
        landsat_match_date_2_sensor <- gsub("L[C|E|T]([0-9]{2}).*", "\\1", basename(landsat_match_date_2_files[1]))
        if(landsat_match_date_2_sensor == "08"){
            landsat_match_date_2_band_list <- list("blue"=2, "green"=3, "red"=4, "nir"=5, "swir1"=6, "swir2"=7)
        }else{
            landsat_match_date_2_band_list <- list("blue"=1, "green"=2, "red"=3, "nir"=4, "swir1"=5, "swir2"=7)
        }
        modis_match_date_2_files <- dir(modis_dir, pattern=paste(".*", strftime(match_date_2, format="%Y_%j"), "_band_[0-9]{1}.*envi$", sep=""), full=T)
        if(length(modis_match_date_2_files) != length(band_names)){
            print(paste("Did not find", length(band_names), "MODIS files for", match_date_2, "in", modis_dir, "|| Aborting"))
            return(NA)
        }
        modis_match_date_2_bands <- as.integer(gsub(".*_band_([0-9]{1})_.*", "\\1", basename(modis_match_date_2_files)))
    }

    # get MODIS file for prediction date
    modis_pred_date_files <- dir(modis_dir, pattern=paste(".*", strftime(pred_date, format="%Y_%j"), "_band_[0-9]{1}.*envi$", sep=""), full=T)
    if(length(modis_pred_date_files) != length(band_names)){
        print(paste("Did not find", length(band_names), "MODIS files for", pred_date, "in", modis_dir, "|| Aborting"))
        return(NA)
    }
    modis_pred_date_bands <- as.integer(gsub(".*_band_([0-9]{1})_.*", "\\1", basename(modis_pred_date_files)))

    ex_r <- raster(landsat_match_date_1_files[1])
    sys_cmds <- c()
    out_files <- c()
    for(band in band_names){
        # determine which files to use for this band
        landsat_file_1 <- landsat_match_date_1_files[which(landsat_match_date_1_bands == as.integer(landsat_match_date_1_band_list[band]))]
        mod_file_1 <- modis_match_date_1_files[which(modis_match_date_1_bands == as.integer(modis_band_list[band]))]
        if(!is.na(match_date_2)){
            landsat_file_2 <- landsat_match_date_2_files[which(landsat_match_date_2_bands == as.integer(landsat_match_date_2_band_list[band]))]
            mod_file_2 <- modis_match_date_2_files[which(modis_match_date_2_bands == as.integer(modis_band_list[band]))]            
        }
        mod_file_pred <- modis_pred_date_files[which(modis_pred_date_bands == as.integer(modis_band_list[band]))]
        out_file <- file.path(out_dir, paste(out_prefix, strftime(pred_date, format="%Y%j"), paste(strftime(Sys.time(), format="%Y%j%H%M"), sep=""), band, sep="_"))
        starfm_file <- file.path(out_dir, paste(basename(out_file), ".txt", sep=""))
        file_con <- file(starfm_file, "w")
        writeLines("STARFM_PARAMETER_START", con=file_con)
        if(!is.na(match_date_2)){
            writeLines("NUM_IN_PAIRS = 2", con=file_con)
            writeLines(paste("IN_PAIR_MODIS_FNAME =", mod_file_1, mod_file_2), con=file_con)
            writeLines(paste("IN_PAIR_LANDSAT_FNAME =", landsat_file_1, landsat_file_2), con=file_con)
        }else{
            writeLines("NUM_IN_PAIRS=1", con=file_con)
            writeLines(paste("IN_PAIR_MODIS_FNAME =", mod_file_1), con=file_con)
            writeLines(paste("IN_PAIR_LANDSAT_FNAME =", landsat_file_1), con=file_con)
        }
        writeLines(paste("IN_PDAY_MODIS_FNAME =", mod_file_pred), con=file_con)
        writeLines(paste("OUT_PDAY_LANDSAT_FNAME =", out_file), con=file_con)
        writeLines(paste("NROWS =", nrow(ex_r)), con=file_con)
        writeLines(paste("NCOLS =", ncol(ex_r)), con=file_con)
        writeLines(paste("RESOLUTION =", res(ex_r)[1]), con=file_con)
        writeLines("SCALE_FACTOR = 10000", con=file_con)
        writeLines("LANDSAT_FILLV = 32767", con=file_con)
        writeLines("LANDSAT_DATA_RANGE = 0, 10000", con=file_con)
        writeLines(paste("LANDSAT_UNCERTAINTY =", landsat_uncertainty), con=file_con)
        writeLines("MODIS_FILLV = 32767", con=file_con)
        writeLines("MODIS_DATA_RANGE = 0, 10000", con=file_con)
        writeLines(paste("MODIS_UNCERTAINTY =", modis_uncertainty), con=file_con)
        if(use_spatial){
            writeLines("USE_SPATIAL_FLAG = 1", con=file_con)
        }else{
            writeLines("USE_SPATIAL_FLAG = 0", con=file_con)
        }
        if(is.na(max_search_dist)){
            writeLines(paste("MAX_SEARCH_DISTANCE =", max(nrow(ex_r), ncol(ex_r))), con=file_con)            
        }else{
            writeLines(paste("MAX_SEARCH_DISTANCE =", max_search_dist), con=file_con)
        }
        writeLines(paste("NUM_SLICE_PURE_TEST =", num_slice), con=file_con)
        writeLines("STARFM_PARAMETER_END", con=file_con)
        close(file_con)

        sys_cmd <- paste(shQuote(normalizePath(starfm_path)), starfm_file)
        if(run_cmd) system(sys_cmd)
        sys_cmds <- c(sys_cmds, sys_cmd)
        out_files <- c(out_files, out_file)
    }
    if(run_cmd){
        return(out_files)
    }else{
        return(sys_cmds)
    }
}

# testing
ard_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard"
mcd43_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43"
ard_vrts <- dir(ard_dir, pattern=".*vrt$", full=T)
which(ard_dates == as.Date("2018-6-26"))
ard_dates <- as.Date(gsub(".*_[0-9]{6}_([0-9]{8})_[0-9]{8}_.*", "\\1", basename(ard_vrts)), format="%Y%m%d")
MakeStarFM(match_date_1=ard_dates[3], match_date_2=ard_dates[62], pred_date=as.Date("2018-6-26"), landsat_dir=ard_dir, modis_dir=mcd43_dir, out_dir=normalizePath("~/Desktop"), out_prefix="starfm_cv_test", run_cmd=TRUE)

# plotting
library(tools)
s_starfm <- stack(sapply(band_names, FUN=function(x) dir(out_dir, pattern=paste("starfm_cv_test.*", x, "$", sep=""), full=T)))
s_actual <- stack(ard_vrts[20])
layout(matrix(1:2, nrow=1))
plotRGB(s_actual, 4, 3, 2, stretch="lin")
plotRGB(s_starfm, 4, 3, 2, stretch="lin")

############################
# Determine possible input pairs
# stupid comment
library(raster)
ard_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_ard"
mcd43_dir <- "/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/central_valley/cv_mcd43"
ard_qas <- dir(ard_dir, pattern="*PIXELQA.*envi$", full=T)
ard_dates <- as.Date(gsub(".*_[0-9]{6}_([0-9]{8})_[0-9]{8}_.*", "\\1", basename(ard_qas)), format="%Y%m%d")
ard_sensors <- as.integer(gsub("L[C|E|T]([0-9]{2}).*", "\\1", basename(ard_qas)))
# we limit possibilities to OLI
ard_qas <- ard_qas[ard_sensors == 8]
ard_dates <- ard_dates[ard_sensors == 8]
# next, compute we'll compute the missing fraction in each
ard_good_fracs <- rep(NA, length(ard_qas))
i <- 1
for(ard_qa in ard_qas){
    print(paste("Doing", basename(ard_qa), "which is", i, "of", length(ard_qas)))
    r <- raster(ard_qa)
    ard_good_fracs[i] <- sum(values(r) == 0) / ncell(r)
    i <- i + 1
}

cv_pair_df <- data.frame(file=ard_qas, date=ard_dates, good_frac=ard_good_fracs)
# This file moved to /Volumes/.../SEAL/KF_four_site_landsat!!!!
# save(cv_pair_df, file="~/Desktop/good_frac_cv.Rdata")
# load("~/Desktop/good_frac_cv.Rdata")

# now, we'll compute the correlation between MODIS images on all possible dates, for each
for(ard_date in ard_dates){
    # find the pred day 
    mod_pred_day <- 
}
