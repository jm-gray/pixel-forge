#!/bin/csh
#$ -pe omp 1
#$ -l h_rt=48:00:00
alias modld 'module load'
module purge
source /project/modislc/src/modis/c5.1/MCD12Q2/bin2hdf_code/prep_environment.sh
time bin2hdf.exe
module purge
source $HOME/.module
