# to submit:
# bsub -q cnr -W 12:00 -n 20 -R "span[ptile=20]" -o /share/jmgray2/KF_output/kf_error.out.%J -e /share/jmgray2/KF_output/kf_error.err.%J R CMD BATCH --vanilla /share/jmgray2/KF_output/KF_error.R


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
cl <- makeCluster(20)
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
system.time(tmp <- parApply(cl, Y[error_sample, ], 1, ProgressiveMissingFraction, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, cdl_process_sds=cdl_process_sds, cdl_types=cdl_types, miss_iter=miss_iter))

# this is the 2nd level of hell, need to preallocate...
total_errors <- do.call(rbind, tmp)

# this should be faster, but not very clean
# in practice, running on a cluster with clean memory, it's not faster at all
# errors_dims <- unlist(lapply(tmp, function(x) dim(x)[1]))
# total_errors <- data.frame(matrix(NA, nrow=sum(errors_dims), ncol=dim(tmp[[1]])[2]))
# names(total_errors) <- names(tmp[[1]])
# last_index <- 1
# i <- 1
# for(N in errors_dims){
#   print(paste(last_index,":",(last_index + N - 1)))
#   total_errors[last_index:(last_index + N - 1), ] <- tmp[[i]]
#   i <- i + 1
#   last_index <- last_index + N
# }


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# save output
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
output_directory <- "/share/jmgray2/KF_output"
output_file <- file.path(output_directory, "kf_error_output.Rdata")
save(total_errors, file=output_file)
