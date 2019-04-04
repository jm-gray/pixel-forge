#!/bin/csh
module load R
R --vanilla < /share/jmgray2/INCA_cluster.R --args -tile $1
