#!/bin/sh
# You have to source this script such that the variables are getting set in the session 
# Set number of threads for other programs
export OMP_NUM_THREADS=1
# Cores on Socket 1 as places
export OMP_PLACES=cores
# Spread threads evenly over cores
export OMP_PROC_BIND=spread
