# NOTE: this script was used to find a good paramerer choice for a smooth
# buildup curve for the object level analysis, some of the curves are more
# deterministic than others

# using ClusterManagers
# using Distributed
const product = Iterators.product

include(joinpath(@__DIR__,"setup.jl"))
datadir = joinpath(@__DIR__,"..","data","count_lengths","run_2018-11-26")

params = load_params(joinpath(datadir,"params.jld2"))
params[!,:pindex] .= 1:size(params,1)
settings = joinpath(@__DIR__,"..","src","settings.toml")
settings = TOML.parsefile(settings)
settings["stimulus"]["repeats"] = 24

settings["bandwidth_ratio"]["window"] = 1.0
settings["bandwidth_ratio"]["delta"] = 0.25

# TODO: get the parameters for the top 10 object models
function only(x)
  @assert length(x) == 1
  first(x)
end

function find_models(ranks)
    select = DataFrame(CSV.read(joinpath(datadir,"model_rankings.csv")) )
    sort!(select,:eratio)
    select = select[(select.t_c_m .> 0) .| (select.t_c_a .> 0),:]

    indices = [
      only(select_params(params,
        t_c_m=select[i,:t_c_m],
        t_c_a=select[i,:t_c_a],Δf=6))
      for i in model_ranks
    ]
    result = params[indices,:]
    result[!,:rank] = model_ranks
    result
end

models = find_models(11:20)

writedir = joinpath(@__DIR__,"..","data","buildup_smooth",string(Date(now())))
isdir(writedir) || mkdir(writedir)

N = 250
runs = collect(product(eachrow(models),1:N))

#@distributed (vcat) for ((name,model),Δ,i) in runs
progress = Progress(prod(size(runs)))
results = mapreduce(vcat,runs) do (model,i)
  with_logger(NullLogger()) do
    results = bistable_model(model,settings,intermediate_results=true)
    len,val = results.percepts.counts
    next!(progress)
    DataFrame(
      length=len,
      rank=model.rank,
      response=val.+1,
      run=i,
    )
  end
end

CSV.write(joinpath(writedir,"build_object_find_smooth_11_20.csv"),results)
