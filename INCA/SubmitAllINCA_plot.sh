#!/bin/bash

# submit all CONUS tiles as INCA tile summary jobs on henry2
SCRATCHDIR="/share/jmgray2/"
SCRIPTNAME="sub_INCA_plot.sh"
SCRATCHROOT="INCA"
RUNTIME="1:00"
declare -a tiles=("h08v04" "h09v04" "h10v04" "h11v04" "h12v04" "h13v04" "h08v05" "h09v05" "h10v05" "h11v05" "h12v05" "h09v06" "h10v06")
# declare -a tiles=("h11v05" "h08v04" "h09v04" "h10v04" "h09v06" "h10v06")
for i in "${tiles[@]}"
do
  bsub_cmd="bsub -q cnr -W $RUNTIME -o ${SCRATCHDIR}${SCRATCHROOT}.$i.out.%J -e $SCRATCHDIR$SCRATCHROOT.$i.err.%J \"csh ${SCRATCHDIR}${SCRIPTNAME} $i\""
  #echo $bsub_cmd
  eval $bsub_cmd
done
