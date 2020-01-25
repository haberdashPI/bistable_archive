# this is a basic test of the key functionality of this project. It is run with
# a minimum working example of the stimulus to veriyf that the code runs
# correctly before running a full simulation on the cluster

using Pkg; Pkg.activate(joinpath(@__DIR__,".."))
include("count_lengths.jl")

# set datadir to the parameter set you want to test
datadir = joinpath(homedir(),"work","dlittle",
                   "bistable_individual_sensitive_W","run_2019-01-02")
# datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-15")

count_lengths(
  datadir=joinpath(datadir,"data"),
  logfile=joinpath(datadir,"logs","test.log"),
  first_index=1,last_index=1,
  sim_repeat=2,
  git_hash="UNKNOWN",
  params=joinpath(datadir,"params.jld2"),
  progressbar=false
)

