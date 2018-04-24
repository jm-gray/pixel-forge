# For Bill Monahan: pixel-scale 0, 5, 25, 50, 75, 95, 100%, 2001-2015, mid-greenup/down, CONUS
library(raster)
library(rgdal)
library(parallel)

# bsub -q cnr -n 12 -W 24:00 -o /rsstu/users/j/jmgray2/SEAL/MonahanQuants.out.%J -e /rsstu/users/j/jmgray2/SEAL/MonahanQuants.err.%J R CMD BATCH --vanilla /rsstu/users/j/jmgray2/SEAL/MonahanQuants.R /rsstu/users/j/jmgray2/SEAL/MonahanQuants.Rlog

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetTileSDS <- function(tile, sds, data_dir, years){
  patt <- paste("(", paste(sds, collapse="|"), ")_", tile, "_(", paste(years, collapse="|"), ").hdf$", sep="")
  return(dir(data_dir, pattern=patt, full=T))
}

GetQuants <- function(x, output_dir="/rsstu/users/j/jmgray2/SEAL/MonahanQuants", data_dir="/gpfs_common/share02/jmgray2/INCA/mcd12q2", years=2001:2015, quants=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)){
  sds <- x[1]
  tile <- x[2]
  in_files <- GetTileSDS(tile, sds, data_dir, years)
  s <- stack(in_files)
  s <- subset(s, seq(1, nlayers(s), by=2))
  NAvalue(s) <- 32767
  doy_offsets <- as.integer(as.Date(paste(years, "-1-1", sep="")) - as.Date("1970-1-1"))
  s <- s - doy_offsets
  v <- values(s)
  system.time(v_quants <- apply(v, 1, function(x) quantile(x, quants, na.rm=T)))
  v_notmissing <- apply(v, 1, function(x) sum(!is.na(x)))
  out_s <- do.call(stack, replicate(length(quants), raster(s, 1)))
  values(out_s) <- t(v_quants)
  out_s_notmissing <- raster(s, 1)
  values(out_s_notmissing) <- v_notmissing
  out_s <- stack(out_s, out_s_notmissing)
  output_file <- file.path(output_dir, paste("Quants_", sds, "_", tile, "_", paste(years[1], years[length(years)], sep="_to_"), ".tif", sep=""))
  writeRaster(out_s, file=output_file)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
data_dir <- "/gpfs_common/share02/jmgray2/INCA/mcd12q2"
sds <- c("MidGreenup", "MidGreendown")
quants <- c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1)
conus_tiles <- c("h08v04", "h09v04", "h10v04", "h11v04", "h12v04", "h13v04", "h08v05", "h09v05", "h10v05", "h11v05", "h12v05", "h08v06", "h09v06", "h10v06")
years <- 2001:2015
NUM_CORES <- 12

apply_df <- expand.grid(sds, conus_tiles)

cl <- makeCluster(NUM_CORES)
clusterExport(cl, c("GetTileSDS", "GetQuants"))
clusterEvalQ(cl, {library(raster); library(rgdal)})
trash <- parApply(cl, apply_df, 1, GetQuants)
