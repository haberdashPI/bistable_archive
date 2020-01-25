using ClusterManagers
using Distributed
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

# if endswith(gethostname(),".cluster")
#     addprocs(SlurmManager(20), partition="CPU", t="24:00:00",
#              cpus_per_task=1)
#     @everywhere include(joinpath(@__DIR__,"..","src","setup.jl"))
# end

models = Dict(
  :object => begin
    p = copy(params[select_params(params,t_c_a=5,t_c_m=5,Δf=6),:])
    p.t_c_σ .= 100.0
    p
  end,
  :central => begin
    p = copy(params[select_params(params,s_c_a=5,s_c_m=5,Δf=6),:])
    p.s_c_σ .= 1.0
    p
  end,
  :peripheral => begin
    p = copy(params[select_params(params,f_c_a=15,f_c_m=130,Δf=6),:])
    p.f_c_σ .= 1.0
    p
  end,
  :combined => begin
    p = copy(params[select_params(params,f_c_a=15,f_c_m=130,Δf=6),:])
    p.f_c_σ .= 1.0
    p.s_c_a .= 5
    p.s_c_m .= 5
    p.s_c_σ .= 1.0
    p.t_c_a .= 5
    p.t_c_m .= 5
    p.t_c_σ .= 1.0
    p
  end
)

writedir = joinpath(@__DIR__,"..","data","buildup",string(Date(now())))
isdir(writedir) || mkdir(writedir)

# run the buildup simulation
N = 1000
deltas = [3,6,12]
runs = collect(product(models,deltas,1:N))

results = @distributed (vcat) for ((name,model),Δ,i) in runs
  # parameter setup
  p = copy(model)
  p[:Δf] = Δ

  with_logger(NullLogger()) do
    results = bistable_model(p,settings,intermediate_results=true)
    len,val = results.percepts.counts
    DataFrame(
      length=len,
      response=val.+1,
      run=i,
      delta_f=Δ,
      level=name
    )
  end
end

CSV.write(joinpath(writedir,"build_results_object_only.csv"),results)
