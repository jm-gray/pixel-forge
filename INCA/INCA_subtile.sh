#!/bin/csh

conda activate /usr/local/usrapps/jmgray2/env_gdal
setenv LD_LIBRARY_PATH /usr/local/usrapps/jmgray2/env_gdal/lib:$LD_LIBRARY_PATH
module load R

R --vanilla < /rsstu/users/j/jmgray2/SEAL/INCA/INCA_cluster.R --args -tile $1
