# This file can be run to generate a new set of parameters to evaluate the
# model on.

# change this parameter to do a different kind of parameter search.
# The kinds of searchs available are listed at the bottom of this file.
kind = :sensitive_W

using Pkg; Pkg.activate(joinpath(@__DIR__,".."))
using FileIO
using DataFrames
using JLD2
using CategoricalArrays
using Unitful
using Unitful: ms, s, Hz, kHz
using Dates
using CorticalSpectralTemporalResponses
using CorticalSpectralTemporalResponses: cycoct

include("parameters.jl")

drop = Base.Iterators.drop

in_ms(x) = Float64.(ustrip.(uconvert.(ms,x)))
in_Hz(x) = Float64.(ustrip.(uconvert.(Hz,x)))

# factoral expansion of parameters
function byparams(params)
  if length(params) == 1
    key,vals = first(params)
    DataFrame(;key => vals)
  else
    key,vals = first(params)
    others = byparams(drop(params,1))

    result = vcat((others for i in 1:length(vals))...)
    result[key] = repeat(vals,inner=nrow(others))
    result
  end
end

unit_table = Dict(r"τ" => in_ms,r"Δt" => in_ms, r"^f$" => in_Hz)

# save all specified parameters under a given label (to be run by `run_all_count_lengths.sh`)
function write_params(label,df)
  dfsize = size(df,1)
  unique!(df)
  @assert size(df,1) == dfsize "All parameter combinations must be unique."

  if occursin("Mycroft",gethostname())
    dir = joinpath(@__DIR__,"..","data",
                   "count_lengths",label,"run_$(Date(now()))")
  elseif occursin("deus1",gethostname())
    dir = joinpath(@__DIR__,"..","data",
                   "count_lengths",label,"run_$(Date(now()))")
  else
    dir = joinpath(homedir(),"work","dlittle","bistable_$label",
                   "run_$(Date(now()))")
  end
  if !isdir(dir)
    mkpath(joinpath(dir,"logs"))
    mkdir(joinpath(dir,"data"))
  end

  open(joinpath(dir,"$(label)_N.txt"),"w") do f
    println(f,"$(size(df,1))")
  end

  save_params(joinpath(dir,"params.jld2"),df)
end

function clean(x)
  d = floor(Int,log10(x))
  round(x,digits=-(d-1))
end

const default_params = Dict(
    :Δt         => [120ms], # SOA of A to B
    :Δf         => [3,6,12], # semitone difference between A and B
    :f          => [500Hz], # frequency of A

    # naming conventions (these differ somewhat from the paper's notation)
    # a = adaptation
    # m = inhibition
    # σ = noise
    # c = magnitude parameter
    # τ = time constant
    # f = peripheral (f)requency representation
    # s = central frequency (s)cale representation
    # t = object (t)racking representation
    # W_m = inhibition matrix related constants
    # W_m_σ = inhibition breadth (Σ_b in the paper)
    # W_m_c = inhibition strength (θ_b in the paper)
    #    _σ_ϕ = inhibition breadth for prior variance (row 1, column 1 of Σ_b in paper)
    #    _σ_N = inhibition breadth for prior mean (row 2, column 2 of Σ_b in paper)
    :f_c_x      => [3.0],    :f_τ_x     => [500ms],
    :f_c_σ      => [0.0],    :f_τ_σ     => [500ms],
    :f_c_a      => [0.0],   :f_τ_a     => [3s],
    :f_c_m      => [0.0],   :f_τ_m     => [350ms],
    :f_W_m_σ    => [5.6],    :f_W_m_c   => [6.0],

    :s_c_x      => [3.0],    :s_τ_x     => [500ms],
    :s_c_σ      => [0.0],    :s_τ_σ     => [500ms],
    :s_c_a      => [0.0],   :s_τ_a     => [3s],
    :s_c_m      => [0.0],   :s_τ_m     => [350ms],
    :s_W_m_σ    => [15.0],   :s_W_m_c   => [6.0],

    :t_c_x      => [3.0],    :t_τ_x     => [500ms],
    :t_c_σ      => [0.0],    :t_τ_σ     => [500ms],
    :t_c_a      => [0.0],  :t_τ_a     => [3s],
    :t_c_m      => [0.0],  :t_τ_m     => [350ms],
    :t_W_m_σ_t  => [7.0],  :t_W_m_σ_ϕ => [7.0],
    :t_W_m_σ_N  => [3.0],  :t_W_m_c    => [6.0]
  )
Params(pairs...) = merge(default_params,Dict(pairs...))


# look at the variations at each level individually
if kind == :individual
  m_vals = a_vals = [0; clean.(10 .^ range(0.7,stop=4,length=8)); 10 .^ 8]

  write_params("individual_levels",vcat(
    byparams(Params(:f_c_σ => [0.2], :f_c_a => a_vals, :f_c_m => m_vals)),
    byparams(Params(:s_c_σ => [0.2], :s_c_a => a_vals, :s_c_m => m_vals)),
    byparams(Params(:t_c_σ => [0.2], :t_c_a => a_vals, :t_c_m => m_vals))))

# this parameter set surveys all variations across the levels simultaneously,
# but for the 6st case only. A second search will cover all three stimuli across
# any parameter sets that generate plausible bistability for 6st.
# (see **TODO NOTEBOOK** for the generation of the second parameter set)
elseif kind == :survey
  m_vals = a_vals = [0; clean.(10 .^ range(0.7,stop=4,length=4))]

  write_params("survey_interact",
    byparams(Params(:Δf => [6],
                    :f_c_σ => [0.2], :f_c_a => a_vals, :f_c_m => m_vals,
                    :s_c_σ => [0.2], :s_c_a => a_vals, :s_c_m => m_vals,
                    :t_c_σ => [0.2], :t_c_a => a_vals, :t_c_m => m_vals)))
  # TODO: also check changes in the shape of the inhibition pattern
  # (I will need to run some preliminary tests and get a working
  # implementation before I can test this)

# this parameter set examines the sensitivity of the results at each individual
# level to changes in several parameters within each analysis level.
# (noise magnitude, inhibition time constant and adaptation time constant)
elseif kind == :sensitive
  m_vals = a_vals = [0; clean.(10 .^ range(0.7,stop=4,length=4))]
  tests = Dict(:c_σ => [0.4,0.8,1.2],
               :τ_m => [150ms,500ms,1s],
               :τ_a => [1s, 5s, 10s])
  alltests = []
  for prefix in ["f_", "s_", "t_"]
    for (suffix,vals) in tests
      for val in vals
        param = Symbol(prefix*string(suffix))
        # @show param
        # @show float(val)
        p = Params(Symbol(prefix*"c_σ") => [0.2],
                   Symbol(prefix*"c_a") => a_vals,
                   Symbol(prefix*"c_m") => m_vals,
                   param => [float(val)])
        push!(alltests,byparams(p))
      end
    end
  end
  write_params("individual_sensitive",vcat(alltests...))
# doing just the noise alone, because a bug in the older version
# of the code above eliminated it.
elseif kind == :sensitive_noise
  m_vals = a_vals = [0; clean.(10 .^ range(0.7,stop=4,length=4))]
  tests = Dict(:c_σ => [0.4,0.8,1.2])
  alltests = []
  for prefix in ["f_", "s_", "t_"]
    for (suffix,vals) in tests
      for val in vals
        param = Symbol(prefix*string(suffix))
        # @show param
        # @show float(val)
        p = Params(Symbol(prefix*"c_σ") => [0.2],
                   Symbol(prefix*"c_a") => a_vals,
                   Symbol(prefix*"c_m") => m_vals,
                   param => [float(val)])
        push!(alltests,byparams(p))
      end
    end
  end
  write_params("individual_sensitive_noise",vcat(alltests...))
# sensitivity to inhibition breadth
elseif kind == :sensitive_W
  m_vals = a_vals = [0; clean.(10 .^ range(0.7,stop=4,length=4))]
  tests = Dict("f_" => Dict(:f_W_m_σ => [2.8, 11.2]),
               "s_" => Dict(:s_W_m_σ => [7.5, 30.0]),
               "t_" => Dict(:t_W_m_σ_t => [3.5, 14.0],
                            :t_W_m_σ_ϕ => [3.5, 14.0],
                            :t_W_m_σ_N => [1.5, 6.0]))

  alltests = []
  for prefix in ["f_", "s_", "t_"]
    for (param,vals) in tests[prefix]
      for val in vals
        p = Params(Symbol(prefix*"c_σ") => [0.2],
                   Symbol(prefix*"c_a") => a_vals,
                   Symbol(prefix*"c_m") => m_vals,
                   param => [float(val)])
        push!(alltests,byparams(p))
      end
    end
  end
  write_params("individual_sensitive_W",vcat(alltests...))
# incomplete test of buildup (not reported in paper)
elseif kind = :buildup
  alltests =
  for τ_factor in range(0.2,stop=6.0,length=20)
    p = Params(:t_c_a => 5, :t_c_m => 5,
              :t_τ_a => 3s * τ_factor,:t_τ_m => 350ms * τ_factor,
              :t_c_σ => 1.2,
              :s_c_σ => 1.2,
              :f_c_σ => 1.2)
    push!(alltests,byparams(p))
  end
end
