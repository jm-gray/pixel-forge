#-------------------------------------------------------------------------------
DownloadMCD43 <- function(x, out_root_dir, overwrite=F){
  # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
  out_dir <- file.path(out_root_dir, "MCD43A4", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
  existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
  if(length(existing_mcd43a4_files) == 0 | overwrite){
    wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A4/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
    # print(wget_cmd)
    # system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
    system(wget_cmd)
    # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
    existing_mcd43a4_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
  }
  out_dir <- file.path(out_root_dir, "MCD43A2", as.character(x[1]), formatC(as.integer(x[2]), width=3, flag=0))
  existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
  if(length(existing_mcd43a2_files) == 0 | overwrite){
    wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[3]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD43A2/', x[1], '/', formatC(as.integer(x[2]), width=3, flag=0), '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
    # print(wget_cmd)
    # system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
    system(wget_cmd)
    # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
    existing_mcd43a2_files <- dir(out_dir, pattern=as.character(x[3]), full=T)
  }
  return(c(existing_mcd43a4_files, existing_mcd43a2_files))
}

#-------------------------------------------------------------------------------
DownloadMCD12Q1 <- function(x, out_root_dir, overwrite=FALSE){
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(argparse)
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
cmd_args <- arg_parser$parse_args()

# tile_to_do <- "h12v04"
start_date <- as.Date("2013-1-1")
end_date <- as.Date("2013-12-31")
daily_dates <- seq.Date(start_date, end_date, by="day")
download_df <- data.frame(year=as.integer(strftime(daily_dates, format="%Y")), doy=as.integer(strftime(daily_dates, format="%j")), tile=cmd_args$tile)

output_dir <- "/rsstu/users/j/jmgray2/SEAL/MODIS_C6"

trash <- apply(download_df, 1, DownloadMCD43, output_dir)
DownloadMCD12Q1(download_df[1, ], output_dir)

# DownloadMCD43 <- function(x, out_root_dir, overwrite=F)
# DownloadMCD12Q1 <- function(x, out_root_dir, overwrite=F)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Rscript C6_Downloader.R -tile h12v04 &
# Current PID: 25765