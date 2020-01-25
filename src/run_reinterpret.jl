using ArgParse
include(joinpath(@__DIR__,"count_lengths.jl"))

parse_settings = ArgParseSettings()
@add_arg_table parse_settings begin
  "files"
    help = "A sequence of data files we wish to reinterpret using new settings."
    required = true
    arg_type = String
    nargs = '+'
  "--logfile", "-l"
    help = "The file to log information to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","data","count_lengths","run.log")
  "--git_hash"
    help = "The git hash of the source code used to run this simulation."
    required = false
    arg_type = String
    default = "DETECT"
  "--folder"
    help = "The subfolder to store reinterpeted results in"
    required = false
    arg_type = String
    default = "reinterpreted"
  "--settings"
    help = "The file specifying all alogirhtm settings for the simulation."
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"settings.toml")
  "--params"
    help = "The file specifying a databse of parameters explored."
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"params.jld2")
end
args = parse_args(parse_settings)

setup_logging(args["logfile"]) do
  git_hash = args["git_hash"]
  @info "Source code hash: "*(git_hash == "DETECT" ? read_git_hash() : git_hash)

  setfile = args["settings"]
  @info("Reinterpretting source masks using settings in $setfile.")
  settings = AuditoryBistabilityLE.read_settings(setfile)

  paramfile = args["params"]
  @info("Loading parameters from $paramfile.")
  params = load_params(paramfile)

  for file in args["files"]
    resultdir = joinpath(dirname(file),args["folder"])
    @info("Reinterpreted file will be stored in $resultdir.")
    isdir(resultdir) || mkdir(resultdir)

    @info("Reinterpretting file $file.")
    jldopen(file,"r") do infile
      jldopen(joinpath(resultdir,basename(file)),"w") do outfile
        for param in keys(infile)
          if occursin(r"param[0-9]{2}",param)
            @info "Loading results for params `$param`."
            for run in keys(infile[param])
              @info "Interpreting run `$run`."
              entry = infile[param][run]
              old_percepts = entry["percepts"]
              @info("Old interpretation had $(length(old_percepts)) percepts.")

              mask = audiospect(entry["mask"],settings)
              lengths, percepts =
                percept_lengths(mask,params[entry["pindex"],:],settings)

              outfile[join((param,run,"lengths"),"/")] = lengths
              outfile[join((param,run,"percepts"),"/")] = percepts
              @info("New interpretation has $(length(percepts)) percepts.")
            end
          end
        end
      end
    end
    @info("DONE")
  end
end
