# to submit:
# bsub -q cnr -W 6:00 -n 16 -R "span[ptile=16]" -o /share/jmgray2/KF_output/kf_error.out.%J -e /share/jmgray2/KF_output/kf_error.err.%J \"csh KF_error.R\


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Preliminaries
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
.libPaths("/home/jmgray2/R/x86_64-unknown-linux-gnu-library/3.1")
library(parallel)
library(dlm)
library(reshape2)
source("/share/jmgray2/KF_output/KF_functions.R")
load("/share/jmgray2/KF_output/nebraska_kf_workspace.Rdata")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Set up cluster
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
cl <- makeCluster(16)
clusterEvalQ(cl, {source("/share/jmgray2/KF_output/KF_functions.R"); library(dlm); library(reshape2)})

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Parse command line options
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# do the error analysis
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# X <- matrix(rnorm(1e5), nrow=1e3)
# system.time(trash <- parApply(cl, X, 1, function(x) max(x, na.rm=T)))
set.seed(42)
num_to_sample <- 1e4
error_sample <- sample(1:dim(Y)[1], num_to_sample)
miss_iter <- 10
# clusterExport(cl, c("ProgressiveMissingFraction"))
# clusterEvalQ(cl, {library(reshape2)})
system.time(total_errors <- parApply(cl, Y[error_sample, ], 1, ProgressiveMissingFraction, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, cdl_process_sds=cdl_process_sds, cdl_types=cdl_types, miss_iter=miss_iter))
total_errors <- do.call(rbind, total_errors)

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# save output
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
output_directory <- "/share/jmgray2/KF_output"
output_file <- file.path(output_directory, "kf_error_output.Rdata")
save(total_errors, file=output_file)
