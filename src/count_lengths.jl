using Logging
using FileIO
using DataFrames
using AuditoryBistabilityLE
using CorticalSpectralTemporalResponses
using Dates
using Printf
using JLD2
using Feather
using InteractiveUtils
using Unitful
using ProgressMeter
using AxisArrays

include("logger.jl")
include("parameters.jl")

struct CountLength
  ratio::Array{Float64}
  bratio::Array{Float64}
  pindex::Int
  created::DateTime
end

function for_results_in(fn,dir;reinterpret=nothing)
  @showprogress for file in readdir(dir)
    if occursin(r"jld2$",file)
      restream = if reinterpret isa Nothing
        nothing
      else
        jldopen(joinpath(dir,reinterpret,file))
      end

      jldopen(joinpath(dir,file),"r") do stream
        for param in keys(stream)
          if occursin(r"param[0-9]{2}",param)
            for run in keys(stream[param])
              entry = stream[param][run]
              if !(restream isa Nothing)
                reentry = Dict("created"  =>   stream[param][run]["created"],
                               "pindex"   =>   stream[param][run]["pindex"],
                               "lengths"  => restream[param][run]["lengths"],
                               "percepts" => restream[param][run]["percepts"])
                fn(reentry)
              else
                fn(entry)
              end # if restream
            end # for run
          end # if occursin
        end # for param
      end # jldopen

      restream isa Nothing || close(restream)
    end
  end
end

count_lengths(args) = count_lengths(;(Symbol(k) => v for (k,v) in args)...)
function count_lengths(;first_index,last_index,
                       datadir=joinpath(data_dir,"count_lengths"),
                       logfile=joinpath(datadir,"sim.log"),
                       params=joinpath(datadir,"params.jld"),
                       dataprefix="results",
                       git_hash="DETECT",
                       sim_repeat=2,
                       settings=joinpath(@__DIR__,"settings.toml"),
                       progressbar=false)
  setup_logging(logfile) do
    @info("Results will be saved to $datadir")
    count_lengths(first_index,last_index,logfile,datadir,dataprefix,
                  params,git_hash,sim_repeat,
                  settings,progressbar)
  end
end

function read_git_hash()
  olddir = pwd()
  cd(@__DIR__)
  hash = read(`git rev-parse HEAD`,String)
  cd(olddir)

  hash
end

function count_lengths(first_index,last_index,logfile,datadir,dataprefix,
                       params,git_hash,sim_repeat,
                       settings,progressbar)
  dir = abspath(datadir)
  isdir(dir) || mkdir(dir)

  @info "Source code hash: "*(git_hash == "DETECT" ? read_git_hash() : git_hash)
  @info "Loading parameters from "*params
  params = load_params(params)
  if first_index > size(params,1)
    error("First parameter index ($first_index) is outside the range of",
          " parameters.")
  end
  last_index = clamp(last_index,first_index,size(params,1))
  indices = first_index:last_index
  @info "Reading parameters for indices $indices"
  @info "Reading settings from file $settings"

  @assert log10(size(params,2)) < 5
  name = @sprintf("%s_params%05d_%05d.jld2",dataprefix,first_index,
                  last_index)
  filename = joinpath(dir,name)

  for i in indices
    @info "Running simulations for parameter $i"
    # determine the number of simulations to run, accounting for any previously
    # run simulations
    # NOTE: assumes only a single process ever accesses this file
    num_repeats = if isfile(filename)
      jldopen(filename,"r") do file
        p = @sprintf("param%05d",i)
        if haskey(file,p)
          max(0,sim_repeat - length(keys(file[p])))
        else
          sim_repeat
        end
      end
    else
      sim_repeat
    end

    @info "$(sim_repeat - num_repeats) simulations run previosuly."
    @info "Running $(num_repeats) more simulations."

    for repeat in 1:num_repeats
      start_time = now()
      params_dict = Dict(k => params[i,k] for k in names(params))
      result = bistable_model(params_dict,settings, progressbar=progressbar)

      jldopen(filename,"a+") do file
        # if it hasn't yet been done, record the times at which the ratios are
        # computed (their sample rate differs, so we need to interpolate them
        # later on).
        if !haskey(file,"btimes_s")
          file["btimes_s"] = ustrip.(uconvert.(s,times(result.percepts.ratio)))
        end

        entry = @sprintf("param%05d",i)
        count = haskey(file,entry) ? length(keys(file[entry])) : 0
        file[@sprintf("param%05d/run%03d/lengths",i,count)] =
          Array(result.percepts.counts[1])
        file[@sprintf("param%05d/run%03d/percepts",i,count)] =
          Array(result.percepts.counts[2])
        file[@sprintf("param%05d/run%03d/mask",i,count)] =
          AuditoryBistabilityLE.compress(result.primary_source)
        file[@sprintf("param%05d/run%03d/pindex",i,count)] = i
        file[@sprintf("param%05d/run%03d/created",i,count)] = start_time

        @info "Completed simulation run $(count+1) for parameter $i."
      end

      @info "Run yielded ~$(length(result.percepts.counts[1])) percepts."
    end
  end
  @info "DONE"
end

