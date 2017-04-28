# cleaning up duplicate Landsat scene downloads, only difference is processing date
# so we delete the older of the two
# data_dir <- "~/Documents/UrmiaLandsat/jmgray2@ncsu.edu-01202017-155716-115"
# data_dir <- "/Volumes/research/fer/jmgray2/NIP/LandsatUrmia/jmgray2@ncsu.edu-01202017-160407-655/"
in_files <- dir(data_dir, pattern=".*.tar.gz$", full=T)
tm_scenes <- gsub(pattern="(.*)-.*", "\\1", basename(in_files))
proc_year <- gsub(pattern=".*-SC([0-9]{4}).*", "\\1", basename(in_files))
tmp <- as.data.frame(table(tm_scenes))
df <- data.frame(infile=in_files, sceneid=tm_scenes, proc_year=proc_year)
df.merged <- merge(df, tmp, by.x="sceneid", by.y="tm_scenes")
old_files <- df.merged$infile[df.merged$Freq > 1 & df.merged$proc_year == 2016]
# delete all the old files
for(old_file in old_files){
  sys_call <- paste("rm", old_file)
  print(sys_call)
  system(sys_call)
}


#######################################################
tm_in_files <- dir("/Volumes/fer$/jmgray2/NIP/LandsatUrmia/jmgray2@ncsu.edu-01202017-160407-655/", pattern=".*.tar.gz$", full=T)
tm_scenes <- gsub(pattern="LT[4|5]([0-9]{6}).*", "\\1", basename(tm_in_files))
tm_dates <- as.Date(gsub(pattern="LT[4|5][0-9]{6}([0-9]{7})-.*", "\\1", basename(tm_in_files)), format="%Y%j")
tm_df <- data.frame(infile=tm_in_files, scene=tm_scenes, date=tm_dates, sensor="TM")

etm_in_files <- dir("/Volumes/fer$/jmgray2/NIP/LandsatUrmia/jmgray2@ncsu.edu-01202017-155716-115/", pattern=".*.tar.gz$", full=T)
etm_scenes <- gsub(pattern="LE7([0-9]{6}).*", "\\1", basename(etm_in_files))
etm_dates <- as.Date(gsub(pattern="LE7[0-9]{6}([0-9]{7})-.*", "\\1", basename(etm_in_files)), format="%Y%j")
etm_df <- data.frame(infile=etm_in_files, scene=etm_scenes, date=etm_dates, sensor="ETM")

oli_in_files <- dir("/Volumes/fer$/jmgray2/NIP/LandsatUrmia/jmgray2@ncsu.edu-01202017-155420-954/", pattern=".*.tar.gz$", full=T)
oli_scenes <- gsub(pattern="LC8([0-9]{6}).*", "\\1", basename(oli_in_files))
oli_dates <- as.Date(gsub(pattern="LC8[0-9]{6}([0-9]{7})-.*", "\\1", basename(oli_in_files)), format="%Y%j")
oli_df <- data.frame(infile=oli_in_files, scene=oli_scenes, date=oli_dates, sensor="OLI")

all_df <- rbind(tm_df, etm_df, oli_df)

all_dates <- c(unique(etm_dates), unique(tm_dates), unique(oli_dates))

plot(all_dates, all_dates, xlim=c(min(all_dates), max(all_dates)), ylim=c(0,1), type="n", xlab="", ylab="", yaxt="n")
cols <- c(brewer.pal(5, "Reds")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Purples")[4])
abline(v=unique(tm_dates), lwd=0.5, col=cols[1])
abline(v=unique(etm_dates), lwd=0.5, col=cols[2])
abline(v=unique(oli_dates), lwd=0.5, col=cols[3])
legend("topleft", legend=c("TM", "ETM+", "OLI"), col=cols, lty=1, lwd=2, pch=NA, bty="o", bg="white")
# abline(v=seq.Date(as.Date("1990-1-1"), as.Date("2000-1-1"), by="year"), lwd=2)
# why the gap in 93-98?


#######################################################
# Preprocessing
# - get new file (check if already processed)
# - unzip files
# - project to common extent
# - cut to basin extent
# - fMask
# - calc SVI (EVI, NDVI, NDWI, TC components, etc.)
# - restack and save
# - goto next file
