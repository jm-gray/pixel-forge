#!/bin/bash
# Example: qsub -V -pe omp 8 -l h_rt=12:00:00 c6_extractor.sh h12v04

echo Submitting tile $1
R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Extractor.R --args -tile $1 -i 2 5 6 7 8 9 10 11 12 13 14 -i_key 2 -prefix midgup midgdown evi2_min evi2_max evi2_area gup_rsquared gup_missing gup_snowcount gdown_rsquared gdown_missing gdown_snowcount -in_dir /projectnb/modislc/users/joshgray/MCD12Q2C6/output -out_dir /projectnb/modislc/users/joshgray/MCD12Q2C6/extracted

