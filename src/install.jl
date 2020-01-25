# run this to get started with a newly cloned repository
using Pkg
root = joinpath(@__DIR__,"..")
Pkg.activate(root)
Pkg.instantiate()

using TOML

# load the config file
config_file = joinpath(@__DIR__,"..","Config.toml")
if isfile(config_file)
  config = TOML.parsefile(config_file)
  if "data" ∉ keys(config)
    error("Config file missing value for `data`.")
  end
  datadir = config["data"]
else
  error("""
Missing a file named `Config.toml` with the following contents:

data = "[folder where bistable data is located]"

This simulation data can be downloaded from this link

https://osf.io/se795/?view_only=a0daa351467f4b84abe3f244d3aaf24e

There are two tar files, which can be unzipped in the same location to create
the final data folder.

""") end

# create data link
link = joinpath(@__DIR__,"..","data")
if !isdir(link)
  symlink(datadir,link)
  @info "The folder `$link` now links to $datadir"
else
  @info "Directory `$link` has already been created."
end
