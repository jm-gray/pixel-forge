# for CLM submission
tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
getVertTile <- function(x) as.integer(strsplit(x,"v")[[1]][2])
vertNums <- sapply(tiles, getVertTile)
sub_tiles <- tiles[vertNums < 7 & vertNums > 1]
out_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/output"
out_prefix <- "c6protoCLMv3_"
for(tile in sub_tiles){
	# check for output tile existence
	if(file.exists(file.path(out_dir, paste(out_prefix, tile, ".Rdata", sep="")))){
		print(paste("Output exists for tile", tile))
		next
	}

	print(paste("Submitting tile: ", tile))
	sys_cmd <- paste("qsub -V -l mem_total=94G -pe omp 16 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_runner.sh", tile)
	# print(sys_cmd)
	system(sys_cmd)
}

# for CLM extraction
tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
getVertTile <- function(x) as.integer(strsplit(x,"v")[[1]][2])
vertNums <- sapply(tiles, getVertTile)
sub_tiles <- tiles[vertNums < 7 & vertNums > 1]
out_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLMextractv2"
out_prefix <- "twentygup_"
for(tile in sub_tiles){
	# check for existence of output files
	existing_files <- Sys.glob(file.path(out_dir, paste(out_prefix, tile, "*.tif", sep="")))
	if(length(existing_files) > 0){
		print(paste("Output exists for tile", tile))
		next
	}

	print(paste("Submitting tile: ", tile))
	sys_cmd <- paste("qsub -V -pe omp 8 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_extractor.sh", tile)
	# print(sys_cmd)
	system(sys_cmd)
}

# check that extraction worked for all files
CheckExistence <- function(x, out_dir){
	files <- Sys.glob(file.path(out_dir, paste(x[1], x[2], x[3], ".tif", sep="*")))
	if(length(files) > 0){
		return(T)
	}else{
		return(F)
	}
	# ifelse(length(files) > 0, return(T), return(F))
}

library(parallel)
cl <- makeCluster(16)
clusterExport(cl, c("CheckExistence"))

tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
getVertTile <- function(x) as.integer(strsplit(x,"v")[[1]][2])
vertNums <- sapply(tiles, getVertTile)
sub_tiles <- tiles[vertNums < 7 & vertNums > 1]
out_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLMextractv2"
years <- 2000:2014
prefixes <- c("twentygup", "fiftygup", "peak", "fiftygdown", "twentygdown", "gupmiss", "gupsnow", "gdownmiss", "gdownsnow", "evi2int", "evi2max", "evi2min", "guprsq", "gdownrsq")
tmp_df <- expand.grid(prefixes, sub_tiles, years, stringsAsFactors=F)

processed_df <- parApply(cl, tmp_df, 1, CheckExistence, out_dir)

missing_df <- tmp_df[!processed_df, ]
names(missing_df) <- c("metric", "tile", "year")

# resubmit the entire tile (not the most efficient on the computer, but in terms of my time)
for(this_tile in unique(missing_df$tile)){
	print(paste("Submitting tile: ", this_tile))
	sys_cmd <- paste("qsub -V -pe omp 8 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_extractor.sh", this_tile)
	# print(sys_cmd)
	system(sys_cmd)

}

#---------------------------------------
# submit all N America tiles
# NOTE: should check for completion...
tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
sub_tiles <- tiles[grep('h(1[2-5])v01|(09|1[0-5])v02|(09|1[0-4])v03|((0[8-9]|1[0-3])v04|(0[8-9]|1[0-2])v05|(0[7-9]|1[0-1])v06)|(0[8-9]|1[0-2])v07', tiles)]

for(tile in sub_tiles){
	print(paste("Submitting tile: ", tile))
	sys_cmd <- paste("qsub -V -l mem_total=64G -pe omp 16 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_runner.sh", tile)
	# print(sys_cmd)
	system(sys_cmd)
}

# submit the extractor scripts if the .Rdata files have been created
# NOTE: should also make sure not to run tiles for which output already exists...
for(tile in sub_tiles){
	if(file.exists(file.path("/projectnb/modislc/users/joshgray/MCD12Q2C6/output", paste("c6protoCLMv3_", tile, ".Rdata", sep="")))){
		print(paste("Submitting tile: ", tile))
		# for alldates
		# sys_cmd <- paste("qsub -V -pe omp 8 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_extractor_alldates.sh", tile)

		# for segqa
		sys_cmd <- paste("qsub -V -pe omp 8 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_extractor_segqa.sh", tile)
		system(sys_cmd)
	}else{
		print(paste(tile, "not processed yet"))
	}
}

#---------------------------------------
# submit all European tiles
# NOTE: should check for completion...
tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
# sub_tiles <- tiles[grep('h(1[2-5])v01|(09|1[0-5])v02|(09|1[0-4])v03|((0[8-9]|1[0-3])v04|(0[8-9]|1[0-2])v05|(0[7-9]|1[0-1])v06)|(0[8-9]|1[0-2])v07', tiles)] # N America
# sub_tiles <- tiles[grep('h(1[0-2])v07|(09|1[0-3])v08|(09|1[0-4])v09|(1[0-4])v10|(1[1-4])v11|(1[1-3])v12|(1[2-3])v13|(1[3-4])v14', tiles)] # S America
sub_tiles <- tiles[grep('h(1[8-9])v01|(1[8-9]|20)v02|(1[7-9]|20)v03|(1[7-9]|20)v04|(1[7-9]|20)v05', tiles)] # Europe

for(tile in sub_tiles){
	print(paste("Submitting tile: ", tile))
	sys_cmd <- paste("qsub -V -l mem_total=64G -pe omp 16 -l h_rt=12:00:00 /projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/c6_runner.sh", tile)
	# print(sys_cmd)
	system(sys_cmd)
}
