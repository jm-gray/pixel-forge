library(raster)

read_int_file <- function(in_file, dsize=1, nbands=1, s_flag=T, bsq_flag=F) {
  ## get the total dimensions of file = nrow*ncol*dsize*nbands
  f_size <- file.info(in_file)$size
  ## open the file
  f <- file(in_file,"rb")
  ## dsize should be an int() 2 bytes
  tot_size <- (f_size/dsize)
  ## can deal with signed or unsigned integers

  temp <- readBin(f,integer(),n=tot_size,endian= "little",size=dsize,signed=s_flag)

  close(f)

  ## if bsq_flag is 1 then we read by row
  ## else we read by column
   ## re-order the temp array into a matrix
   if(bsq_flag==1) { byr_flag = FALSE } else { byr_flag = TRUE }
      temp <- matrix(temp, ncol=nbands, byrow=byr_flag)

  return(temp)
}

pheno_data_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLMextractv2"
lc_data_dir <- "/projectnb/modislc/data/mcd12_out/lc_out/c5_bin/c5.1_classifications/bug_fix_c5.1_081214/mcd12y3/500m_out/"
metrics <- c("fiftygup", "peak", "fiftygdown", "evi2int", "evi2max", "evi2min")
# tiles <- c("h12v04", "h11v05", "h11v04", "h12v03")
tiles <- c("h11v03", "h12v03", "h10v04", "h11v04", "h12v04", "h10v05", "h11v05")
num_pixels_per_tile <- 1e6
years <- 2001:2014
set.seed(42)

# first_time <- T

# initialize output matrix
sample_mat <- matrix(NA, nrow=length(tiles) * length(years) * num_pixels_per_tile, ncol=length(metrics) + 4)
i <- 1

for(tile in tiles){
  #DEBUG
  print(paste("Working on tile", tile))

  # Get tile's stable 2004-2006 LC
  tile_lc_file <-  dir(lc_data_dir, pattern=paste("2004-2006.igbp1.*", tile, sep=""), full=T)
  lc_v <- read_int_file(tile_lc_file)

  # Get LC proporations for stratified random sample
  lc_breaks <- c(seq(-0.5, 16.5, by=1), 256)
  h_lc <- hist(lc_v, breaks=lc_breaks, plot=F)
  lc_pixel_nums <- round(h_lc$counts[2:17] / sum(h_lc$counts[2:17]) * num_pixels_per_tile)
  sample_lc <- unlist(lapply(1:16, function(x) rep(x, lc_pixel_nums[x])))
  pixel_ids <- 1:5760000
  lc_ids <- lapply(1:16, function(x) pixel_ids[lc_v == x])
  sample_ids <- unlist(lapply(1:length(lc_pixel_nums), function(x) sample(lc_ids[[x]], lc_pixel_nums[x])))

  # check for too many samples due to rounding
  if(length(sample_ids) > num_pixels_per_tile){
    num_elim <- length(sample_ids) - num_pixels_per_tile
    elim_ids <- sample(1:length(sample_ids), num_elim)
    sample_ids <- sample_ids[-1 * elim_ids]
    lc_ids <- lc_ids[-1 * elim_ids]
    sample_lc <- sample_lc[-1 * elim_ids]
  }

  # loop through each year, create a stack of phenology data and retrieve metrics
  for(year in years){
    #DEBUG
    print(paste("Working on year", year))

    in_files <- sapply(metrics, function(x) file.path(pheno_data_dir, paste(x, "_", tile, "_", year, ".tif", sep="")))
    s <- stack(in_files)
    doy_offset <- as.numeric(as.Date(paste(year, "-1-1", sep=""), format="%Y-%j") - as.Date("1970-1-1"))
    s_adj <- c(rep(doy_offset, 3), rep(0, 3))
    s <- s - s_adj
    s_v <- s[sample_ids]

    row_start <- (num_pixels_per_tile * (i - 1) + 1)
    row_end <- row_start + length(sample_ids) - 1
    sample_mat[row_start:row_end,] <- cbind(rep(tile, length(sample_ids)), sample_ids, sample_lc, rep(year, length(sample_ids)), s_v)
    # sample_mat[(num_pixels_per_tile * (i - 1) + 1):(num_pixels_per_tile * i),] <- cbind(rep(tile, length(sample_ids)), sample_ids, sample_lc, rep(year, length(sample_ids)), s_v)
    i <- i + 1

    # if(first_time){
    #   # sample_mat <- cbind(sample_lc, rep(year, length(sample_ids)), s_v)
    #   sample_mat <- cbind(rep(tile, length(sample_ids)), sample_ids, sample_lc, rep(year, length(sample_ids)), s_v)
    #   first_time <- F
    # }else{
    #   # sample_mat <- rbind(sample_mat, cbind(sample_lc, rep(year, length(sample_ids)), s_v))
    #   sample_mat <- rbind(sample_mat, cbind(rep(tile, length(sample_ids)), sample_ids, sample_lc, rep(year, length(sample_ids)), s_v))
    # }
  }
}

# sample_mat <- na.omit(sample_mat)
sample_mat <- as.data.frame(sample_mat)
for(i in 2:10) sample_mat[,i] <- as.numeric(as.character(sample_mat[,i]))
save(sample_mat, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/PCA/seven_tile_sample.Rdata")

# # perform the PCA
# system.time(s_pca <- prcomp(na.omit(sample_mat[, c(-1, -2)]), center=T, scale.=T))
# system.time(s_pca <- prcomp(sample_mat[, c(-1, -2, -3, -4)], center=T, scale.=T)
s_pca <- prcomp(sample_mat[, c(5, 6, 7)], center=T, scale.=T)
summary(s_pca)
s_pca$rotation

PCARaster <- function(pca, tile, years, metrics, pca_cols=1:6, num_components=3, data_dir="/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLMextractv2"){
  first_time <- T
  year_counter <- 1
  for(year in years){
    in_files <- sapply(metrics, function(x) file.path(data_dir, paste(x, "_", tile, "_", year, ".tif", sep="")))
    s <- stack(in_files)
    doy_offset <- as.numeric(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1"))
    s <- s - c(rep(doy_offset, 3), rep(0, 3))
    s_v <- values(s)
    s_v_pca <- predict(pca, newdata=s_v[, pca_cols])
    tmp_r <- raster(s, 1)

    if(first_time){
      out_stack <- lapply(1:num_components, function(x) stack(tmp_r))
      first_time <- F
      for(i in 1:num_components) out_stack[[i]] <- setValues(out_stack[[i]], s_v_pca[, i])
    }else{
      for(i in 1:num_components){
        tmp_pc_r <- setValues(tmp_r, s_v_pca[, i])
        out_stack[[i]] <- stack(out_stack[[i]], tmp_pc_r)
      }
    }
    year_counter <- year_counter + 1
  }
  return(out_stack)
}

metrics <- c("fiftygup", "peak", "fiftygdown", "evi2int", "evi2max", "evi2min")
years <- 2003:2014
h12v04_pca_3comp <- PCARaster(s_pca, tile="h12v04", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h11v05_pca_3comp <- PCARaster(s_pca, tile="h11v05", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h11v04_pca_3comp <- PCARaster(s_pca, tile="h11v04", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h12v05_pca_3comp <- PCARaster(s_pca, tile="h12v05", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h13v04_pca_3comp <- PCARaster(s_pca, tile="h13v04", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h10v05_pca_3comp <- PCARaster(s_pca, tile="h10v05", years=years, metrics=metrics, pca_cols=1:3, num_components=3)
h10v04_pca_3comp <- PCARaster(s_pca, tile="h10v04", years=years, metrics=metrics, pca_cols=1:3, num_components=3)

#====================================================
# trend estimation (see c6_trends.R)
cl <- makeCluster(16)
clusterEvalQ(cl, {library(mblm)})
clusterExport(cl, c("CalcTheilSen"))

pca_result_s <- h12v04_pca_3comp
# pca_result_s <- h11v05_pca_3comp

pca1_v <- values(pca_result_s[[1]])
ts_result_pca1 <- parApply(cl, pca1_v, 1, CalcTheilSen, years)
slope_raster_pca1 <- setValues(tmp_r, ts_result_pca1[1,])
pvalue_raster_pca1 <- setValues(tmp_r, ts_result_pca1[2,])
plot(slope_raster_pca1)

pca2_v <- values(pca_result_s[[2]])
ts_result_pca2 <- parApply(cl, pca2_v, 1, CalcTheilSen, years)
slope_raster_pca2 <- setValues(tmp_r, ts_result_pca2[1,])
pvalue_raster_pca2 <- setValues(tmp_r, ts_result_pca2[2,])
plot(slope_raster_pca2)

pca3_v <- values(pca_result_s[[2]])
ts_result_pca3 <- parApply(cl, pca3_v, 1, CalcTheilSen, years)
slope_raster_pca3 <- setValues(tmp_r, ts_result_pca3[1,])
pvalue_raster_pca3 <- setValues(tmp_r, ts_result_pca3[2,])
plot(slope_raster_pca3)

slope_raster_pca1[pvalue_raster_pca1 > 0.05] <- 0
slope_raster_pca2[pvalue_raster_pca2 > 0.05] <- 0
slope_raster_pca3[pvalue_raster_pca3 > 0.05] <- 0

slope_raster_pca1 <- abs(slope_raster_pca1)
slope_raster_pca2 <- abs(slope_raster_pca2)
slope_raster_pca3 <- abs(slope_raster_pca3)

total_slope <- slope_raster_pca1 + slope_raster_pca2 + slope_raster_pca3


####################################
# get the means, merge, and make a plot
h13v04_pc1_mean <- stackApply(h13v04_pca_3comp[[1]],rep(1,12),mean)
h13v04_pc2_mean <- stackApply(h13v04_pca_3comp[[2]],rep(1,12),mean)
h13v04_pc3_mean <- stackApply(h13v04_pca_3comp[[3]],rep(1,12),mean)
h13v04_pca_mean <- stack(h13v04_pc1_mean, h13v04_pc2_mean, h13v04_pc3_mean)

h12v04_pc1_mean <- stackApply(h12v04_pca_3comp[[1]],rep(1,12),mean)
h12v04_pc2_mean <- stackApply(h12v04_pca_3comp[[2]],rep(1,12),mean)
h12v04_pc3_mean <- stackApply(h12v04_pca_3comp[[3]],rep(1,12),mean)
h12v04_pca_mean <- stack(h12v04_pc1_mean, h12v04_pc2_mean, h12v04_pc3_mean)

h11v05_pc1_mean <- stackApply(h11v05_pca_3comp[[1]],rep(1,12),mean)
h11v05_pc2_mean <- stackApply(h11v05_pca_3comp[[2]],rep(1,12),mean)
h11v05_pc3_mean <- stackApply(h11v05_pca_3comp[[3]],rep(1,12),mean)
h11v05_pca_mean <- stack(h11v05_pc1_mean, h11v05_pc2_mean, h11v05_pc3_mean)

h11v04_pc1_mean <- stackApply(h11v04_pca_3comp[[1]],rep(1,12),mean)
h11v04_pc2_mean <- stackApply(h11v04_pca_3comp[[2]],rep(1,12),mean)
h11v04_pc3_mean <- stackApply(h11v04_pca_3comp[[3]],rep(1,12),mean)
h11v04_pca_mean <- stack(h11v04_pc1_mean, h11v04_pc2_mean, h11v04_pc3_mean)

h12v05_pc1_mean <- stackApply(h12v05_pca_3comp[[1]],rep(1,12),mean)
h12v05_pc2_mean <- stackApply(h12v05_pca_3comp[[2]],rep(1,12),mean)
h12v05_pc3_mean <- stackApply(h12v05_pca_3comp[[3]],rep(1,12),mean)
h12v05_pca_mean <- stack(h12v05_pc1_mean, h12v05_pc2_mean, h12v05_pc3_mean)

h10v05_pc1_mean <- stackApply(h10v05_pca_3comp[[1]],rep(1,12),mean)
h10v05_pc2_mean <- stackApply(h10v05_pca_3comp[[2]],rep(1,12),mean)
h10v05_pc3_mean <- stackApply(h10v05_pca_3comp[[3]],rep(1,12),mean)
h10v05_pca_mean <- stack(h10v05_pc1_mean, h10v05_pc2_mean, h10v05_pc3_mean)

h10v04_pc1_mean <- stackApply(h10v04_pca_3comp[[1]],rep(1,12),mean)
h10v04_pc2_mean <- stackApply(h10v04_pca_3comp[[2]],rep(1,12),mean)
h10v04_pc3_mean <- stackApply(h10v04_pca_3comp[[3]],rep(1,12),mean)
h10v04_pca_mean <- stack(h10v04_pc1_mean, h10v04_pc2_mean, h10v04_pc3_mean)

tmp_merge <- merge(h12v05_pca_mean, h12v04_pca_mean, h11v05_pca_mean, h11v04_pca_mean, h13v04_pca_mean, h10v04_pca_mean, h10v05_pca_mean)
plotRGB(tmp_merge, stretch="lin")
