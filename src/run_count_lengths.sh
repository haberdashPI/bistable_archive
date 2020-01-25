#!/bin/sh
#SBATCH
#SBATCH --job-name=bistable
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --requeue
#SBATCH --cpus-per-task=8

# module load julia
julia  -e 'using Pkg; Pkg.activate("projects/bistable")' \
       -e 'include("projects/bistable/src/run_count_lengths.jl")' \
       -O3 --banner=no $@

