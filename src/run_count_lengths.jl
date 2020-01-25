using ArgParse
include(joinpath(@__DIR__,"count_lengths.jl"))

parse_settings = ArgParseSettings()
@add_arg_table parse_settings begin
  "first_index"
    help = "First parameter index to test."
    required = false
    arg_type = Int
    default = 1
  "last_index"
    help = "Last parameter index to test."
    required = false
    arg_type = Int
    default = 0
  "--git_hash"
    help = "The git hash of the source code used to run this simulation."
    required = false
    arg_type = String
    default = "DETECT"
  "--params"
    help = "The file specifying a databse of parameters to explore."
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"params.jld2")
  "--sim_repeat", "-r"
    help = "How many times to repeat simulation for each parameter."
    required = false
    arg_type = Int
    default = 1
  "--datadir", "-d"
    help = "The directory to save results to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","data","count_lengths")
  "--dataprefix"
    help = "The prefix appended to all result files"
    required = false
    arg_type = String
    default = "results"
  "--logfile", "-l"
    help = "The file to log information to"
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"..","..","data","count_lengths","run.log")
  "--settings"
    help = "The file specifying all alogirhtm settings for the simulation."
    required = false
    arg_type = String
    default = joinpath(@__DIR__,"settings.toml")
end
args = parse_args(parse_settings)

count_lengths(args)
