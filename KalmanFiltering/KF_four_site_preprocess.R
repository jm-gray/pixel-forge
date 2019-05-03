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
# Define the subsets
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
plot(cv_modis_aea, lty=2)
plot(cv_landsat_aea, add=T, lty=1)

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
