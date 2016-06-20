#!/bin/bash
# Example: qsub -V -pe omp 8 -l h_rt=12:00:00 -l mem_total=98G c6_pheno_extractor.sh h09v06

echo Submitting tile $1
R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_NumCycles.R --args -tile $1 -i 3



