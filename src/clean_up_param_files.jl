include(joinpath(@__DIR__,"count_lengths.jl"))
# datadir = joinpath(@__DIR__,"..","data","count_lengths")
datadir = joinpath(homedir(),"work","dlittle","bistable_individual_levels")
for dir in readdir(datadir)
  if occursin(r"run_[0-9-]+",dir) && isdir(joinpath(datadir,dir))
    @info "Updating params in $dir."

    files = readdir(joinpath(datadir,dir))
    param_file = joinpath(datadir,dir,
                          files[findfirst(x -> occursin("params",x),files)])
    @info "Loading params $param_file..."
    params = load_params(param_file)

    oldfile = param_file*".old"
    mv(param_file,oldfile)
    @info "Saved old file to $oldfile."

    newfile = joinpath(datadir,dir,"params.jld2")
    save_params(newfile,params)
    @info "New parameters saved in $newfile"
  end
end
