library(raster)
library(parallel)
library(RColorBrewer)

# directories, constants, and parameters
c5_dir <- "/projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/"
c5_sub_dirs <- c("Greenup", "Maturity", "Senescence", "Dormancy")
c6_dir <- "/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6"
c6_sub_dirs <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
lc_dir <- "/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp"
years <- 2001:2014
lc_year <- 2012
samples_per_tile <- 1e4
sangram_tiles <- c("h12v03", "h12v04", "h20v07", "h20v08", "h23v02", "h27v05")
my_tiles <- c("h11v04", "h11v02", "h12v09", "h25v06", "h27v05")
tiles <- c(sangram_tiles, my_tiles)

# c5 file paths: /projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/2010/001/Greenup/PheGre.A2010001.h12v04.hdf
# c6 file paths: /projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6/Greenup_h12v04_2010
# lc file paths: /projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp/igbp.h12v04.2010.bip

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetStratCellSample <- function(tile, year=2012, num_samples=1e4, exclude_lcs=c(17, 255), lc_range=c(1, 17), lc_dir="/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I12/igbp", rand_seed=42){
  set.seed(rand_seed)
  in_file <- dir(lc_dir, pattern=paste(tile, ".*", year, ".bip$", sep=""), full=T)
  lc_v <- values(raster(in_file))
  lc_tb <- table(lc_v)
  lc_inds <- !as.integer(names(lc_tb)) %in% exclude_lcs
  samples_per_lc <- round(num_samples * lc_tb[lc_inds] / sum(lc_tb[lc_inds]))
  sample_df <- data.frame(lc=as.integer(names(lc_tb)[lc_inds]), lc_samples=as.integer(samples_per_lc))
  pixel_ids <- 1:length(lc_v)
  lc_cell_samples <- apply(sample_df, 1, function(x) sample(pixel_ids[lc_v == x[1]], x[2]))
  names(lc_cell_samples) <-  as.integer(names(lc_tb)[lc_inds])
  return(lc_cell_samples)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetC5Data <- function(tile, year, cell_samples, c5_dir="/projectnb/modislc/data/mcd12_out/phen_out/c5_hdf1/", c5_sub_dirs=c("Greenup", "Maturity", "Senescence", "Dormancy")){
  # gets the C5 MCD12Q2 data for a particular tile-year
  # returns C6 metrics (estimating 50% GUP/GDOWN and peak), cast as UNIX-epoch time
  in_files <- lapply(c5_sub_dirs, function(x) dir(file.path(c5_dir, year, "001", x), pattern=paste(tile, ".*.hdf$", sep=""), full=T))
  s <- stack(in_files)
  v <- values(s)
  doy_offset <- as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("2000-1-1") )
  doy_add <- as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1") )
  v <- v - doy_offset + doy_add - 1 # convert to UNIX epoch (days since 1970-1-1)
  # estimate additional C6 SDS
  midgreenup <- (v[, 1:2] + v[, 3:4]) / 2
  peak <- (v[, 3:4] + v[, 5:6]) / 2
  midgreendown <- (v[, 5:6] + v[, 7:8]) / 2
  v <- cbind(v, midgreenup, peak, midgreendown)
  v <- v[, c(1:2, 9:10, 3:4, 11:12, 5:6, 13:14, 7:8)] # rearrange in phenometric order
  # retrieve the values for the specified cells
  pheno_samples <- do.call(rbind, lapply(cell_samples, function(x) v[x, ]))
  # package up as a data.frame with  and return
  lc_ids <- rep(as.integer(names(cell_samples)), sapply(cell_samples, length))
  pheno_df <- data.frame(tile=rep(tile, length(lc_ids)), year=rep(year, length(lc_ids)), lc=lc_ids, cell_num=unlist(cell_samples))
  pheno_names <- c("c5greenup", "c5midgreenup", "c5maturity", "c5peak", "c5senescence", "c5midgreendown", "c5dormancy")
  pheno_df <- cbind(pheno_df, pheno_samples)
  names(pheno_df)[5:18] <- paste(rep(pheno_names, each=2), c(1, 2), sep="_")
  return(pheno_df)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetC6Data <- function(tile, year, cell_samples, c6_dir="/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6", c6_sub_dirs=c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")){
  in_files <- lapply(c6_sub_dirs, function(x) dir(file.path(c6_dir, x), pattern=paste(tile, ".*", year, "$", sep=""), full=T))
  s <- stack(in_files)
  NAvalue(s) <- 32767
  v <- values(s)
  # retrieve the values for the specified cells
  pheno_samples <- as.data.frame(do.call(rbind, lapply(cell_samples, function(x) v[x, ])))
  pheno_names <- c("c6greenup", "c6midgreenup", "c6maturity", "c6peak", "c6senescence", "c6midgreendown", "c6dormancy")
  names(pheno_samples) <- paste(rep(pheno_names, each=2), c(1, 2), sep="_")
  return(pheno_samples)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetTileYearData <- function(tile, year){
  cell_samples <- GetStratCellSample(tile=tile, year=year)
  c5_df <- GetC5Data(tile=tile, year=year, cell_samples=cell_samples)
  c6_df <- GetC6Data(tile=tile, year=year, cell_samples=cell_samples)
  tile_year_df <- cbind(c5_df, c6_df)
  return(tile_year_df)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# run over all tile year combinations
all_tile_years <- expand.grid(tiles, years)
cl <- makeCluster(16)
clusterEvalQ(cl, {library(raster)})
clusterExport(cl, c("GetTileYearData", "GetStratCellSample", "GetC6Data", "GetC5Data"))
all_pheno_data_list <- apply(all_tile_years, 1, function(x) GetTileYearData(tile=x[1], year=x[2]))
c5_c6_pheno_data <- do.call(rbind, all_pheno_data_list)
save(c5_c6_pheno_data, file="/projectnb/modislc/users/joshgray/C6_Diagnostics/c5_c6_intercomparison.Rdata")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plotting
MakePhenometricDensityPlot <- function(subdf, phenometric_name="Greenup"){
  # df should be two column: 1st column the C5 data, second the C6
  c5_doys <- as.integer(strftime(as.Date(subdf[,1], origin="1970-1-1"), format="%j"))
  c6_doys <- as.integer(strftime(as.Date(subdf[,2], origin="1970-1-1"), format="%j"))
  lims <- range(c(c5_doys, c6_doys), na.rm=T)
  pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  smoothScatter(c5_doys, c6_doys, xlim=lims, ylim=lims, colramp=pal, nrpoints=0, xlab="C5 DOY", ylab="C6 DOY")
  title(phenometric_name)
  abline(a=0, b=1, lwd=2, lty=2, col=1)
}

MakePhenometricScatterByLCPlot <- function(df, cols, phenometric_name="Greenup"){
  # df should be 3 column: 1st column is LC ID, 2nd column the C5 data, 3rd the C6
  plot_cols <- cols[df[,1]]
  c5_doys <- as.integer(strftime(as.Date(df[,2], origin="1970-1-1"), format="%j"))
  c6_doys <- as.integer(strftime(as.Date(df[,3], origin="1970-1-1"), format="%j"))
  lims <- range(c(c5_doys, c6_doys))
  plot(c5_doys, c6_doys, col=plot_cols, xlab="C5 DOY", ylab="C6 DOY", pch=16, cex=0.25)
  # pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  # smoothScatter(c5_doys, c6_doys, xlim=lims, ylim=lims, colramp=pal, nrpoints=0, xlab="C5", ylab="C6")
  title(phenometric_name)
}

# igbp colors and codes
igbp_colors <- c("#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97", "#8FBB8F", "#BB8F8F", "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766", "#4682B2", "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD", "#000000")
igbp_names <- c("water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs", "open shrubs", "woody savannas", "savannas", "grasslands", "perm wetlands", "croplands", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified")
igbp_codes <- c(17, 1:16, 255)
MakeTransp <- function(x, alpha=255/2) rgb(t(col2rgb(x)), alpha=alpha, max=255)
igbp_col_df <- data.frame(igbp_codes, igbp_names, igbp_colors, igbp_t_cols=MakeTransp(igbp_colors))[order(igbp_codes),]
col_vec <- rep(NA, max(igbp_col_df$igbp_codes))
col_vec[igbp_col_df$igbp_codes] <- as.character(igbp_col_df$igbp_t_cols)

# make density plots of overall agreement for cycle 1 for each phenometric
c5_c6_pheno_data <- c5_c6_pheno_data[c5_c6_pheno_data$lc != 255, ] # didn't exclude unclassified the first time through!
layout(matrix(1:8, nrow=2, byrow=T))
pheno_names <- c("greenup", "midgreenup", "maturity", "peak", "senescence", "midgreendown", "dormancy")
c5_cols <- seq(5, 18, by=2)
c6_cols <- seq(19, 32, by=2)
i <- 1
for(phenometric in pheno_names){
  MakePhenometricDensityPlot(c5_c6_pheno_data[, c(c5_cols[i], c6_cols[i])], phenometric_name=pheno_names[i])
  i <- i + 1
}

# make scatter plots of overall agreement for cycle 1 for each phenometric, colorized by LC
layout(matrix(1:8, nrow=2, byrow=T))
pheno_names <- c("greenup", "midgreenup", "maturity", "peak", "senescence", "midgreendown", "dormancy")
c5_cols <- seq(5, 18, by=2)
c6_cols <- seq(19, 32, by=2)
i <- 1
for(phenometric in pheno_names){
  MakePhenometricScatterByLCPlot(c5_c6_pheno_data[, c(3, c5_cols[i], c6_cols[i])], phenometric_name=pheno_names[i], cols=col_vec)
  i <- i + 1
}

# MakePhenometricScatterByLCPlot(c5_c6_pheno_data[, c(3, 5, 19)], phenometric_name=pheno_names[i], cols=col_vec)

phenometric_missing_fraction <- apply(c5_c6_pheno_data[, 5:32], 2, sum(is.na(x)) / length(x))
