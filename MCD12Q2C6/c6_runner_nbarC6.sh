#!/bin/bash
# Example: qsub -V -l mem_total=128G -pe omp 16 -l h_rt=12:00:00 c6_runner_nbarC6.sh h12v04

echo Submitting tile $1
R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Cluster.R --args -tile h12v04 -data_dir /projectnb/modislc/users/dsm/nbars_c6 -C6input -gup 0.2 0.5 1 -gdown 0.2 0.5 -min_seg_amplitude 0.15 -start_date 2013-1-1 -spline_spar 0 -out_prefix pheno_c6nbar_ -out_dir /projectnb/modislc/users/joshgray/MCD12Q2C6/output_nbar_C6

