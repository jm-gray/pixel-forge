SCRATCHDIR="/rsstu/users/j/jmgray2/SEAL/INCA/"
SCRIPTNAME="INCA_subtile.sh"
SCRATCHROOT="INCA"

declare -a tiles=("h08v04" "h09v04" "h10v04" "h11v04" "h12v04" "h13v04" "h08v05" "h09v05" "h10v05" "h11v05" "h12v05" "h08v06" "h09v06" "h10v06")
for i in "${tiles[@]}"
do
  bsub_cmd="bsub -q cnr -W 4:00 -n 16 -R \"oc span[ptile=16]\" -o ${SCRATCHDIR}${SCRATCHROOT}.$i.out.%J -e $SCRATCHDIR$SCRATCHROOT.$i.err.%J \"csh ${SCRATCHDIR}${SCRIPTNAME} $i\""
  echo $bsub_cmd
#   eval $bsub_cmd
done

