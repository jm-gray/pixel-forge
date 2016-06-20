#!/bin/bash
# Example: qsub -V -l mem_total=128G -pe omp 16 -l h_rt=12:00:00 c6_runner.sh h12v04

echo Submitting tile $1
R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Cluster.R --args -tile $1 -gup 0.1 0.5 1 -gdown 0.1 0.5 -min_seg_amplitude 0.1 -out_prefix c6proto_


