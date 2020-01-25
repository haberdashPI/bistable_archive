import Distributions: TDist, Normal, logpdf, pdf

abstract type Stats{T} end

update(stats::Stats{T},x::AbstractVector{T}) where T =
  update!(copy(stats),x)

include(joinpath(@__DIR__,"priors","beta.jl"))
include(joinpath(@__DIR__,"priors","ridge_norm.jl"))
include(joinpath(@__DIR__,"priors","grouping.jl"))
