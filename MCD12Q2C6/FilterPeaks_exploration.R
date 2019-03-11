library(raster)
library(rgdal)
library(parallel)
library(RColorBrewer)

source("~/Documents/pixel-forge/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
DownloadMCD12Q2C6 <- function(x, out_root_dir, overwrite=F){
  print(paste("Downloading:", as.character(x[2]), x[1]))
  if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
  wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[2]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/1467/MCD12Q2/', x[1], '/001/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_root_dir, "/", x[1], sep="")
  system(wget_cmd)
}

cl
tiles <- c("h08v05", "h08v07", "h11v04", "h11v09", "h12v04", "h12v09", "h12v12", "h13v12", "h18v02", "h19v02", "h25v03", "h25v04", "h25v06", "h26v05", "h27v05", "h31v11")
all_tiles <- read.table("~/Google Drive/DataSets/gltiles.txt")
new_tiles <- as.character(all_tiles$V1[!(all_tiles$V1 %in% tiles)])
# download_df <- expand.grid(2001:2016, tiles)
download_df <- expand.grid(2001:2016, new_tiles)
data_dir <- "/Volumes/users/j/jmgray2/SEAL/MCD12Q2C6"
apply(download_df, 1, DownloadMCD12Q2C6, out_root_dir=data_dir)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
DownloadTile <- function(x, out_root_dir, overwrite=F){
  # x is c(year, tile, product)
  print(paste("Downloading:", as.character(x[2]), x[1], x[3]))
  out_file <- file.path(out_root_dir, x[3])
  wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[2]), '*tar.gz" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/', x[3], '/', x[1], '/001/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_root_dir, "/", x[3], sep="")
  # wget -e robots=off -m -np -nd -R -r -l1 -A "*h08v05*tar.gz" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD12I4/2012/001/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P /Volumes/rsstu/users/j/jmgray2/SEAL/MCD12I4
  # if(file.exists(out_file) & !overwrite){
  #   print(paste(out_file, "exists, skipping"))
  # }else{
  #   system(wget_cmd)
  # }
  system(wget_cmd)
  # print(wget_cmd)
}

download_df <- expand.grid(2012:2015, "h08v05", c("MCD12I2", "MCD12I3", "MCD12I4"))
data_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/"
apply(download_df, 1, DownloadTile, out_root_dir=data_dir)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetCellTimeSeries <- function(cell, r, evi2_files, resid_files, snow_files, year_to_do=2009){
  the_line <- rowFromCell(r, cell)
  the_col <- colFromCell(r, cell)
  tmp <- Get3YearDataChunk_tmp(evi2_files, resid_files, snow_files, year_to_do, the_line, 1)
  return(tmp[the_col,])
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get a collection of cell time series on which to test the algorithm
data_dir="~/Desktop/MCD12Q2_analysis"
splined_nbar_data_dir <- file.path(data_dir, "MCD12I2")
splined_qa_data_dir <- file.path(data_dir, "MCD12I3")
splined_resids_data_dir <- file.path(data_dir, "MCD12I4")

tiles <- gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", dir(splined_nbar_data_dir, pattern="bip"))
years <- gsub(".*A([0-9]{4}).*", "\\1", dir(splined_nbar_data_dir, pattern="bip"))
year_to_do <- 2009
c6_dates <- as.Date(paste(rep((year_to_do - 1):(year_to_do + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")
# tiles_i2 <- gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", dir(splined_nbar_data_dir, pattern="bip"))
# years_i2 <- gsub(".*A([0-9]{4}).*", "\\1", dir(splined_nbar_data_dir, pattern="bip"))
# tiles_i3 <- gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", dir(splined_qa_data_dir, pattern="bip"))
# years_i3 <- gsub(".*A([0-9]{4}).*", "\\1", dir(splined_qa_data_dir, pattern="bip"))
# tiles_i4 <- gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", dir(splined_resids_data_dir, pattern="bip"))
# years_i4 <- gsub(".*A([0-9]{4}).*", "\\1", dir(splined_resids_data_dir, pattern="bip"))

GetTileSamples <- function(tile_to_do, nsamples=100, year_to_do=2009, lw_mask_dir="/Volumes/research/fer/jmgray2/MODIS/LWMASK500/"){
  metric <- "Greenup"
  r <- raster(dir(lw_mask_dir, pattern=paste(".*", tile_to_do, ".*bin$", sep=""), full=T))
  # r_v1 <- raster(dir(phen_dir_v1, pattern=paste("^", metric, ".*", tile_to_do, ".*", year_to_do, "$", sep=""), full=T))
  # NAvalue(r_v1) <- 32767
  tmp_df <- data.frame(id=1:ncell(r), v=values(r))
  tmp_df <- subset(tmp_df, tmp_df$v == 1)
  cells_to_do <- sample(tmp_df$id, nsamples)

  patt <- paste("MCD12", ".*", tile_to_do, ".*.bip$", sep="")
  evi2_files <- dir(splined_nbar_data_dir, pattern=patt, full=T)
  resid_files <- dir(splined_resids_data_dir, pattern=patt, full=T)
  snow_files <- dir(splined_qa_data_dir, pattern=patt, full=T)

  tile_ts <- t(sapply(cells_to_do, GetCellTimeSeries, r=r, evi2_files=evi2_files, resid_files=resid_files, snow_files=snow_files))
  tile_ts <- cbind(data.frame(tile=tile_to_do, cell=cells_to_do), tile_ts)
  return(tile_ts)
}

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("GetTileSamples", "GetCellTimeSeries", "splined_nbar_data_dir", "splined_qa_data_dir", "splined_resids_data_dir", "Get3YearDataChunk_tmp", "ReadSplinedNBAR"))
clusterEvalQ(cl, {library(raster); library(rgdal)})
system.time(test_ts <- do.call(rbind, parLapply(cl, unique(tiles), GetTileSamples, nsamples=1000)))
test_ts[test_ts == 32767] <- NA
save(test_ts, file="/Volumes/research/fer/jmgray2/MODIS/test_ts.Rdata")

# do the phenology ala v1 and updated rel_amp
pars_v1 <- DefaultPhenoParameters()
pars_v1$min_seg_amplitude <- 0.1

pars_rel_amp <- DefaultPhenoParameters()
pars_rel_amp$min_seg_amplitude <- 0.1
pars_rel_amp$rel_amp_frac <- 0.35
pars_rel_amp$rel_peak_frac <- 0.7

pars_rel_amp_180 <- DefaultPhenoParameters()
pars_rel_amp_180$min_seg_amplitude <- 0.1
pars_rel_amp_180$rel_amp_frac <- 0.35
pars_rel_amp_180$rel_peak_frac <- 0.7
pars_rel_amp_180$max_increase_length <- 180
pars_rel_amp_180$max_decrease_length <- 180

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(raster); library(rgdal); source("~/Documents/pixel-forge/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")})

pheno_result_v1 <- parApply(cl, test_ts[, 3:3287], 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pars_v1, pheno_period_start=as.Date("2009-1-1"), pheno_period_end=as.Date("2009-12-31"))
pheno_result_rel_amp <- parApply(cl, test_ts[, 3:3287], 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pars_rel_amp, pheno_period_start=as.Date("2009-1-1"), pheno_period_end=as.Date("2009-12-31"))
pheno_result_rel_amp_180 <- parApply(cl, test_ts[, 3:3287], 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pars_rel_amp_180, pheno_period_start=as.Date("2009-1-1"), pheno_period_end=as.Date("2009-12-31"))

pheno_result_v1_m <- matrix(pheno_result_v1, ncol=25, byrow=T)
pheno_result_rel_amp_m <- matrix(pheno_result_rel_amp, ncol=25, byrow=T)
pheno_result_rel_amp_180_m <- matrix(pheno_result_rel_amp_180, ncol=25, byrow=T)

pheno_results <- cbind(test_ts[,1:2], pheno_result_v1_m, pheno_result_rel_amp_m, pheno_result_rel_amp_180_m)
pheno_names <- c("num_cyles", paste(rep(c("evi_area", "evi_amp", "evi_min", "ogi", "midgup", "mat", "peak", "sen", "midgdown", "dor", "oqa", "dqa"), 2), rep(1:2, each=12), sep="_"))
names(pheno_results) <- c("tile", "cell", paste(pheno_names, "v1", sep="_"), paste(pheno_names, "rel_amp", sep="_"), paste(pheno_names, "rel_amp_180", sep="_"))
pheno_results[pheno_results == 32767] <- NA

PlotDiffHist <- function(x, metric="ogi", comp1="v1", comp2="rel_amp", which_cycle=1){
  comp1_col <- which(names(x) == paste(metric, which_cycle, comp1, sep="_"))
  comp2_col <- which(names(x) == paste(metric, which_cycle, comp2, sep="_"))
  hist(x[, comp1_col] - x[, comp2_col])
}


which_series <- 4
trash_v1 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars_v1, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
trash_v1_rel_amp <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars_rel_amp, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
trash_v1_rel_amp_180 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars_rel_amp_180, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)

layout(matrix(1:3, nrow=3))
par(mar=c(0, 4, 0, 2), oma=c(4, 0, 1, 1))
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(trash_v1[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[1], lwd=1.5)
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(trash_v1_rel_amp[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(trash_v1_rel_amp_180[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[3])

# map of select tiles
continents <- shapefile("~/Downloads/Cont/continent_sin.shp")
mod_grid <- shapefile("~/Downloads/modis_grid (2)/modis_sinusoidal_grid_world.shp")
oceans <- shapefile("~/Downloads/Cont/OceansSin.shp")
mod_grid$tile <- paste("h", formatC(as.integer(mod_grid$h), width=2, flag="0"), "v", formatC(as.integer(mod_grid$v), width=2, flag="0"), sep="")
plot(oceans, col=brewer.pal(5, "Blues")[2], border=NA)
plot(continents, col=grey(0.4), border=NA, add=T)
plot(mod_grid, border=grey(0.7), add=T)
plot(mod_grid[mod_grid$tile %in% as.character(unique(test_ts$tile)), ], border=2, add=T, lwd=1.5)

v_edges <- seq(extent(mod_grid)[1], extent(mod_grid)[2], len=37)
h_edges <- seq(extent(mod_grid)[3], extent(mod_grid)[4], len=19)
# text(x=v_edges[1:36] + diff(v_edges) / 2, y=max(h_edges), labels=1:36, pos=3, col=1)
text(x=v_edges[1:36] + diff(v_edges) / 2, y=min(h_edges), labels=0:35, pos=1, col=1)
text(x=min(v_edges), y=h_edges[1:18] + diff(h_edges) / 2, labels=17:0, pos=2, col=1)
# abline(v=v_edges[1:36] + diff(v_edges) / 2)
# abline(h=h_edges[1:18] + diff(h_edges) / 2)





# h11v09: 75, 65, 51, 18
# h12v09: 75
# h12v12: 89, 35, 23, 68 <- a strong indicator that we shouldn't use the relative peak?
# h18v02: 25
# h25v06: 42, 26, 85 (71 is great dbl cropping)
# h31v11: 16, 6, 88
# h25v04: 46


#TESTING
mycols <- c(brewer.pal(5, "Reds")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Greens")[4], brewer.pal(5, "Purples")[4])
pars_v1 <- DefaultPhenoParameters()
pars_v1$min_seg_amplitude <- 0.1

pars_rel_amp <- DefaultPhenoParameters()
pars_rel_amp$min_seg_amplitude <- 0.1
pars_rel_amp$rel_amp_frac <- 0.35
pars_rel_amp$rel_peak_frac <- 0.7

pars_rel_amp_180 <- DefaultPhenoParameters()
pars_rel_amp_180$min_seg_amplitude <- 0.1
pars_rel_amp_180$rel_amp_frac <- 0.35
# pars_rel_amp_180$rel_peak_frac <- 0.7
pars_rel_amp_180$max_increase_length <- 180
pars_rel_amp_180$max_decrease_length <- 180

pars_rel_amp_no_peak <- DefaultPhenoParameters()
pars_rel_amp_no_peak$min_seg_amplitude <- 0.1
pars_rel_amp_no_peak$rel_amp_frac <- 0.35


n <- sample(1:1e3, 1)
x <- as.numeric(subset(test_ts, tile == "h11v09")[n, 3:ncol(test_ts)])

trash_v1 <- AnnualPhenologyC6(as.numeric(x), c6_dates, pars_v1, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
trash_v1_rel_amp <- AnnualPhenologyC6(as.numeric(x), c6_dates, pars_rel_amp, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
trash_v1_rel_amp_180 <- AnnualPhenologyC6(as.numeric(x), c6_dates, pars_rel_amp_180, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
trash_v1_rel_amp_no_peak <- AnnualPhenologyC6(as.numeric(x), c6_dates, pars_rel_amp_no_peak, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)

layout(matrix(1:4, nrow=4))
par(mar=c(0, 4, 0, 2), oma=c(4, 0, 1, 1))
PlotSeries(x, c6_dates)
abline(v=as.Date(trash_v1[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[1], lwd=1.5)
PlotSeries(x, c6_dates)
abline(v=as.Date(trash_v1_rel_amp[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])
PlotSeries(x, c6_dates)
abline(v=as.Date(trash_v1_rel_amp_180[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[3])
PlotSeries(x, c6_dates)
abline(v=as.Date(trash_v1_rel_amp_no_peak[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[4])


#@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$
#@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$
# All old below here
#@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$
#@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
lc_dir <- "~/Desktop/MCD12Q2_analysis/LC/"
data_dir="~/Desktop/MCD12Q2_analysis"
splined_nbar_data_dir <- file.path(data_dir, "MCD12I2")
splined_qa_data_dir <- file.path(data_dir, "MCD12I3")
splined_resids_data_dir <- file.path(data_dir, "MCD12I4")
phen_dir_v0 <- file.path(data_dir, "MCD12I6")
phen_dir_v1 <- file.path(data_dir, "MCD12I6_v1")

tile_to_do <- "h08v05"
year_to_do <- 2009
metric <- "Greenup"
c6_dates <- as.Date(paste(rep((year_to_do - 1):(year_to_do + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Select some places where dormancy is delayed in v1
r_v0 <- raster(dir(phen_dir_v0, pattern=paste("^", metric, ".*", tile_to_do, ".*", year_to_do, "$", sep=""), full=T))
r_v1 <- raster(dir(phen_dir_v1, pattern=paste("^", metric, ".*", tile_to_do, ".*", year_to_do, "$", sep=""), full=T))
NAvalue(r_v0) <- 32767
NAvalue(r_v1) <- 32767

r_diff <- r_v0 - r_v1
v_diff <- values(r_diff)
# not_zero <- which(v_diff != 0)
# cells_to_do <- not_zero[c(1, 3, 9, 7, 12, 14, 15, 17, 26)] # hand selected
cells_to_do <- sample(which(v_diff <0 & !is.na(v_diff)), 5)
v_diff[cells_to_do]
# cell_lines <- rowFromCell(r_v0, cells_to_do)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# get the phenology data for v0
pheno_in_files <- dir(phen_dir_v0, pattern=paste("(Greenup|MidGreenup|Maturity|Peak|Senescence|MidGreendown|Dormancy).*", tile_to_do, "_[0-9]{4}$", sep=""), full=T)
pheno_years <- as.integer(gsub(".*([0-9]{4})$", "\\1", basename(pheno_in_files)))
pheno_in_files <- pheno_in_files[pheno_years %in% c(year_to_do - 1, year_to_do, year_to_do + 1)]
reorder_inds <- c(4:6, 13:15, 7:9, 16:18, 19:21, 10:12, 1:3)
pheno_s_v0 <- stack(pheno_in_files[reorder_inds])
NAvalue(pheno_s_v0) <- 32767
pheno_v_v0 <- values(pheno_s_v0)

# get the phenology data for v1
pheno_in_files <- dir(phen_dir_v1, pattern=paste("(Greenup|MidGreenup|Maturity|Peak|Senescence|MidGreendown|Dormancy).*", tile_to_do, "_[0-9]{4}$", sep=""), full=T)
pheno_years <- as.integer(gsub(".*([0-9]{4})$", "\\1", basename(pheno_in_files)))
pheno_in_files <- pheno_in_files[pheno_years %in% c(year_to_do - 1, year_to_do, year_to_do + 1)]
reorder_inds <- c(4:6, 13:15, 7:9, 16:18, 19:21, 10:12, 1:3)
pheno_s_v1 <- stack(pheno_in_files[reorder_inds])
NAvalue(pheno_s_v1) <- 32767
pheno_v_v1 <- values(pheno_s_v1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# retrieve the splined input data
patt <- paste("MCD12", ".*", tile_to_do, ".*.bip$", sep="")
evi2_files <- dir(splined_nbar_data_dir, pattern=patt, full=T)
resid_files <- dir(splined_resids_data_dir, pattern=patt, full=T)
snow_files <- dir(splined_qa_data_dir, pattern=patt, full=T)

tile_ts <- t(sapply(cells_to_do, GetCellTimeSeries, r=r_v0, evi2_files=evi2_files, resid_files=resid_files, snow_files=snow_files))

mycols <- c(brewer.pal(5, "Reds")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Greens")[4])
layout(matrix(1:nrow(tile_ts), ncol=1, byrow=T))
par(mar=c(1, 4, 1, 1))
for(i in 1:nrow(tile_ts)){
  # PlotSeries(tile_ts[i,], c6_dates, igbp_names)
  PlotSeries(tile_ts[i,], c6_dates, plot_legend=F)
  abline(v=as.Date(pheno_v_v0[cells_to_do[i],], origin=as.Date("1970-1-1")), lty=rep(1:2, 21), col=mycols[1], lwd=1.5)
  abline(v=as.Date(pheno_v_v1[cells_to_do[i],], origin=as.Date("1970-1-1")), lty=rep(1:2, 21), col=mycols[2])
  legend("topleft", legend=c("v0", "v1"), lty=1, col=mycols[1:2])
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
which_series <- 5
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.15
pars$rel_amp_frac <- 0.35
pars$rel_peak_frac <- 0.67
# pars$rel_peak_frac <- 0.7
pheno_results_v0 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.1
pheno_results_v1 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
layer_names <- c("num_cycles", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "overall_qa_cycle1", "detailed_qa_cycle1",  "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "overall_qa_cycle2", "detailed_qa_cycle2")

ind_pattern <- c(rep(NA, 4), 3, 9, 15, 21, 27, 33, 39, rep(NA, 5), 4, 10, 16, 22, 28, 34, 40, NA, NA)
data.frame(name=layer_names, v0_mine=pheno_results_v0, v0_giss=pheno_v_v0[cells_to_do[which_series],][ind_pattern], v1_mine=pheno_results_v1, v1_giss=pheno_v_v1[cells_to_do[which_series],][ind_pattern])

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# test combo of min_seg_amplitude and rel_amp
which_series <- 5
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.1
pheno_only_min_seg <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars$rel_amp_frac <- 0.35
pheno_with_rel_amp <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars$rel_peak_frac <- 0.7
pheno_with_rel_amp_and_rel_peak <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)

# data.frame(layer_names, no_filt=pheno_no_filter, filt=pheno_filtered)
layout(matrix(1:3, nrow=3))
par(mar=c(0, 4, 0, 2), oma=c(4, 0, 1, 1))
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_only_min_seg[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[1], lwd=1.5)
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_with_rel_amp[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_with_rel_amp_and_rel_peak[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Testing the peak-to-peak distance problem

x[x==pheno_pars$nbar_NA_value] <- NA
data_length <- length(x) / 3
evi <- x[1:data_length]
evi <- evi / pheno_pars$nbar_scale_factor
resids <- x[(data_length + 1):(2 * data_length)]
resids <- resids / pheno_pars$nbar_scale_factor
snowflags <- x[(2 * data_length + 1):(3 * data_length)]

which_series <- 3
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.1
pheno_no_filter <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars$min_peak_to_peak_distance=90
pheno_filtered <- AnnualPhenologyC6_filter(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
# data.frame(layer_names, no_filt=pheno_no_filter, filt=pheno_filtered)
layout(matrix(1:2, nrow=2))
par(mar=c(0, 4, 0, 2), oma=c(4, 0, 1, 1))
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_no_filter[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[1], lwd=1.5)
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_filtered[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])

FilterPeaks <- function(peaks, min_peak_to_peak_distance=NA){
  # this version takes advantage of the fact that the peaks are in magnitude order, so no need to check EVI
  peaks <- rev(peaks)
  i <- 1
  while(i < length(peaks)){
    if(!is.na(peaks[i])){
      peaks_too_close <- which(abs(peaks[i] - peaks) < min_peak_to_peak_distance)
      peaks_too_close <- peaks_too_close[peaks_too_close != i]
      peaks[peaks_too_close] <- NA
      print(peaks)
    }
    i <- i + 1
  }
  peaks <- peaks[!is.na(peaks)]
  if(length(peaks) == 0) return(NA)
  return(peaks)
}


# FilterPeaks <- function(x, potential_peaks, min_peak_to_peak_distance=NA){
#
#   potential_peaks <- sort(potential_peaks)
# 	# check that potential_peaks is not NA
# 	if(all(is.na(potential_peaks))) return(NA)
#
#   # return unmodified peaks if the min_peak_to_peak_distance is NA or 0
#   if(is.na(min_peak_to_peak_distance) | min_peak_to_peak_distance == 0) return(potential_peaks)
#
# 	# we loop through all the potential peaks, only moving on when we have satisfied
# 	# the minimum distance between peaks requirement
# 	i <- 1
#
# 	while(i < length(potential_peaks)){
#     # print(paste(i, potential_peaks[i]))
# 		# find the distance to next peak
# 		peak_to_peak_dist <- potential_peaks[i + 1] - potential_peaks[i]
#
# 		if(peak_to_peak_dist < min_peak_to_peak_distance){
# 			# eliminate the smaller of the two potential peaks
# 			if(x[potential_peaks[i + 1]] >= x[potential_peaks[i]]){
# 				potential_peaks <- potential_peaks[-1 * i]
# 			}else{
# 				potential_peaks <- potential_peaks[-1 * (i + 1)]
# 			}
# 		}else{
# 			# distance is fine, move on
# 			i <- i + 1
# 		}
# 		# if we've eliminated the last peak, return NA
# 		if(length(potential_peaks) == 0) return(NA)
# 	}
# 	return(potential_peaks)
# }

AnnualPhenologyC6_filter <- function(x, dates, pheno_pars, pheno_period_start, pheno_period_end, plot=F){
	# processes MCD12Q2 phenology for a single pixel's time series
	# Up to two cycles are returned (but total # of peaks in period is too)
	# if there were >2 segs in the period, then the two with the highest magnitude
	# peaks are returned, but in chronological order

	# get a template return value
	annual_pheno_metrics <- PhenoReturnValue()

	# split the data up: first third are smoothed nbar-evi2, second third are residuals, and last third are snowflags
	# also scale, and fill for NA
	x[x==pheno_pars$nbar_NA_value] <- NA
	data_length <- length(x) / 3
	evi <- x[1:data_length]
	evi <- evi / pheno_pars$nbar_scale_factor
	resids <- x[(data_length + 1):(2 * data_length)]
	resids <- resids / pheno_pars$nbar_scale_factor
	snowflags <- x[(2 * data_length + 1):(3 * data_length)]

	# find valid peaks in the time series
	# valid_peaks <- try(FilterPeaks(evi, FindPotentialPeaks_minmax(evi), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	valid_peaks <- try(FindPeaks(evi), silent=T)
  valid_peaks <- FilterPeaks(valid_peaks, min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance)
	# if(inherits(valid_peaks, 'try-error') | is.na(valid_peaks)){
	if(inherits(valid_peaks, 'try-error') | (length(valid_peaks) == 0)){
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# find full segments
	full_segs <- try(GetSegs(valid_peaks, evi, pheno_pars), silent=T)
	if(inherits(full_segs, 'try-error') | is.null(full_segs)){
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# eliminate segs that end before, or start after the period of interest
	seg_overlaps <- unlist(lapply(full_segs, SegOverlapsPeriod, dates, pheno_period_start, pheno_period_end))
	if(all(!seg_overlaps)){
		# no segs within period of interest
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}
	full_segs <- full_segs[seg_overlaps]

	# get the segment metrics
	seg_metrics <- lapply(full_segs, SegMet, x=x, dates=dates, pheno_pars=pheno_pars)

	# limit to segments where the peak is in the period of interest
	is_in_period <- unlist(lapply(seg_metrics, PeakInPeriod, start_date=pheno_period_start, end_date=pheno_period_end))
	segs_in_period <- seg_metrics[is_in_period]

	# assign the seg metric values to the output return value
	num_cycles <- length(segs_in_period) # total number of valid peaks within the period of interest
	if(num_cycles == 0){
		# no segments in period of interest, return default annual pheno metrics
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}else{
		# so, there is at least one valid segment.
		segs_in_period <- rev(segs_in_period) # reverse order so highest magnitude segs are first
		segs_in_period <- segs_in_period[order(unlist(segs_in_period)[which(names(unlist(segs_in_period)) == "peak")])] # put into chronological order

		# set first cycle return metrics
		cycle1_seg_met <- segs_in_period[[1]]
		annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle1_seg_met, cycle=1, num_cycles=num_cycles, fill_code=0)
		if(num_cycles > 1){
			# set the second cycle
			cycle2_seg_met <- segs_in_period[[2]]
			annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle2_seg_met, cycle=2)
		}
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}
}
