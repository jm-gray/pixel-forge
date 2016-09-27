#!/bin/csh
source ./prep_environment.sh
make -f makefile_bin2hdf
make -f makefile_bin2hdf_phe_temp
make -f makefile_bin2hdf_phe
module purge
source $HOME/.module
