#!/bin/csh

source /usr/local/apps/R/R-312.csh
R --vanilla < /share/jmgray2/INCA_plot.R --args -tile $1
