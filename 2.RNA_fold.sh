#!/bin/bash


#SBATCH -J print_st #Slurm job name

# Set the maximum runtime, uncomment if you need it
##SBATCH -t 48:00:00 #Maximum runtime of 48 hours

# Choose partition (queue), for example, partition "standard"
##SBATCH -p standard
#SBATCH -p general
##SBATCH -p ssci
##SBATCH -p himem

# Use 2 nodes and 48 cores
#SBATCH -N 1 -n 1

RNAfold --noPS < TLR10-15-reference.fa > TLR10-15-reference.fa.dotstring
