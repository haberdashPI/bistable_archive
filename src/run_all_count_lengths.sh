#!/bin/sh

# runs a defined search of parameter space: see params_for_count_lengths.jl
# to generate a new search of parameter space

# command line arguments (only the first two are required)
result_dir=$1 # where to store results, and where to read the parameter set from, e.g. data/count_lengths/run_2018-11-09/
label=$2 # the label to use for this simulation, (ala `write_params`, as defined in params_for_count_lengths.jl)
K=${3:-10} # the number of parameter sets to run per node in the cluster
repeat=${4:-20} # the number of repetitions of the simulation to run
default_N=`cat ${result_dir}/${label}_N.txt` # the number of parameter configurations (generated by params_for_count_lengths.jl)
N=${5:-$default_N}

# set this to the directory where this file is stored
proj_dir="projects/bistable/src"

# what iteratio to start on? (change to resume failed runs)
S=1

cd ${proj_dir}
GIT_HASH=`git rev-parse HEAD`
cd

# this just echos the commands, once you verify that it's right, pipe it to sh
for i in `seq $S $K $N`; do
  echo "sbatch ${proj_dir}/run_count_lengths.sh $i \
    $((i+K-1)) -r ${repeat} --git_hash $GIT_HASH \
    --params ${result_dir}/params.jld2 \
    --settings ${proj_dir}/settings.toml \
    -d ${result_dir}/data/ \
    -l ${result_dir}/logs/result_${i}.log"
done

