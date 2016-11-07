#!/bin/csh

source /usr/local/apps/R/R-312.csh
#R CMD BATCH --vanilla /share/jmgray2/INCA_cluster.R --args -tile $1
R --vanilla < /share/jmgray2/INCA_cluster.R --args -tile $1

