#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Preprocesses MCD43A4/2 data: reprojects & resamples to match an example raster
# cut to cutline, calculate EVI2. Applies to R,G,B,NIR and QA, & QAsnow
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(raster)
library(rgdal)
library(tools)
library(parallel)

GetSDSName <- function(mcd43a4_file_path, band){
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
}

GetSDSNameQA <- function(mcd43a2_file_path){
	return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality", sep = ""))
}

GetSDSNameSnow <- function(mcd43a2_file_path){
	return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
}

GetBandBRDFQualities <- function(x){
	# represent the binary conversion of x as a 32 x N matrix of values
	M <- matrix(as.integer(intToBits(t(x))), ncol=32, byrow=T)
	M <- M[, 32:1] # reverse column order

	# get the individual bit-words and convert to decimal
	QA_Fill <- M[1]
	NOT_USED <- M[2:4]
	B7QA <- sum(M[5:8] * t(c(8,4,2,1)))
	B6QA <- sum(M[9:12] * t(c(8,4,2,1)))
	B5QA <- sum(M[13:16] * t(c(8,4,2,1)))
	B4QA <- sum(M[17:20] * t(c(8,4,2,1)))
	B3QA <- sum(M[21:24] * t(c(8,4,2,1)))
	B2QA <- sum(M[25:28] * t(c(8,4,2,1)))
	B1QA <- sum(M[29:32] * t(c(8,4,2,1)))

	return(c(QA_Fill, B1QA, B2QA, B3QA, B4QA, B5QA, B6QA, B7QA))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ModisOveralapPreprocessNBAR <- function(modis_file, cutline_file, example_raster_file, suffix="_modis_overlap", out_dir=NA){
  # TM Band 1 is blue, MODIS equivalent is Band 3
  # TM Band 2 is green, MODIS equivalent is Band 4
  # TM Band 3 is red, MODIS equivalent is Band 1
  # TM Band 4 is NIR, MODIS equivalent is Band 2
  tmp_r <- raster(example_raster_file)
  bands <- 1:4
  band_names <- c("red", "nir", "blue", "green")
  for(i in 1:length(bands)){
    out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(modis_file)), "_", band_names[i], suffix, ".tif", sep=""))
    gdal_cmd <- paste("gdalwarp -overwrite -dstnodata -9999 -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(tmp_r), "'", sep=""), "-te", paste(extent(tmp_r)[c(1,3,2,4)], collapse=" "), "-ts", ncol(tmp_r), nrow(tmp_r), "-crop_to_cutline -cutline", cutline_file, GetSDSName(modis_file, bands[i]), out_file)
    system(gdal_cmd)
  }

  # calc evi2
  red_file <- file.path(out_dir, paste(file_path_sans_ext(basename(modis_file)), "_", band_names[1], suffix, ".tif", sep=""))
  nir_file <- file.path(out_dir, paste(file_path_sans_ext(basename(modis_file)), "_", band_names[2], suffix, ".tif", sep=""))
  evi2_out_file <- gsub("red", "evi2", red_file)
  gdal_cmd <- paste("gdal_calc.py -A ", nir_file, " -B ", red_file, " --calc='2.5 * (((0.0001 * A) - (0.0001 * B)) / ((0.0001 * A) + (2.4 * (B * 0.0001)) + 1.0)) * 10000' --NoDataValue=-9999 --overwrite --outfile=", evi2_out_file, sep="")
  system(gdal_cmd)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ModisOveralapPreprocessQA <- function(modis_file, cutline_file, example_raster_file, suffix="_modis_overlap", out_dir=NA){
  tmp_r <- raster(example_raster_file)

  # process the BRDF Albedo Band Quality
  out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(modis_file)), "_", "band_quality", suffix, ".tif", sep=""))
  gdal_cmd <- paste("gdalwarp -overwrite -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(tmp_r), "'", sep=""), "-te", paste(extent(tmp_r)[c(1,3,2,4)], collapse=" "), "-ts", ncol(tmp_r), nrow(tmp_r), "-crop_to_cutline -cutline", cutline_file, GetSDSNameQA(modis_file), out_file)
  system(gdal_cmd)

  # process the BRDF Snow Albedo
  out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(modis_file)), "_", "snow", suffix, ".tif", sep=""))
  gdal_cmd <- paste("gdalwarp -overwrite -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(tmp_r), "'", sep=""), "-te", paste(extent(tmp_r)[c(1,3,2,4)], collapse=" "), "-ts", ncol(tmp_r), nrow(tmp_r), "-crop_to_cutline -cutline", cutline_file, GetSDSNameSnow(modis_file), out_file)
  system(gdal_cmd)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# run script over all directories to create Landsat KF input
# get all the ETM+ and TM data directories
example_raster_file <- "/projectnb/modislc/users/joshgray/DL_Landsat/LT50290312010361EDC00/LT50290312010361EDC00_sr_band1_landsat_overlap.tif"
cutline_file <- "/projectnb/modislc/users/joshgray/DL_Landsat/landsat_overlap.shp"
tile <- "h10v04"

cl <- makeCluster(16)
clusterExport(cl, c("ModisOveralapPreprocessNBAR", "GetSDSName", "GetSDSNameQA", "GetSDSNameSnow", "ModisOveralapPreprocessQA"))
clusterEvalQ(cl, {library(tools); library(raster); library(rgdal)})

# process the NBAR data
nbar_in_files <- c()
for(year in 2007:2010){
  data_dir <- paste("/projectnb/modislc/data/mcd12_in/c5/mcd43a4/", year, sep="")
  nbar_in_files <- c(nbar_in_files, dir(data_dir, pattern=paste(".*", tile, ".*hdf$", sep=""), full=T, rec=T))
}

# apply the Stack and Subset function to each directory
trash <- parLapply(cl, nbar_in_files, ModisOveralapPreprocessNBAR, example_raster_file=example_raster_file, cutline_file=cutline_file, out_dir="/projectnb/modislc/users/joshgray/DL_Landsat/MODISOVERLAP")

# process the NBAR QA data
qa_in_files <- c()
for(year in 2007:2010){
  data_dir <- paste("/projectnb/modislc/data/mcd12_in/c5/mcd43a2/", year, sep="")
  qa_in_files <- c(qa_in_files, dir(data_dir, pattern=paste(".*", tile, ".*hdf$", sep=""), full=T, rec=T))
}

# apply the Stack and Subset function to each directory
trash <- parLapply(cl, qa_in_files, ModisOveralapPreprocessQA, example_raster_file=example_raster_file, cutline_file=cutline_file, out_dir="/projectnb/modislc/users/joshgray/DL_Landsat/MODISOVERLAP")
