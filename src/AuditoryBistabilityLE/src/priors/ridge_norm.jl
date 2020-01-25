export ridgenorm, rdist

function rdist(;scale,freq)
  ((as,af),(bs,bf)) -> (as - bs)^2 / (scale^2) + (af - bf)^2 / (freq^2)
end

logdet_(x,n) = logdet(x)
logdet_(x::UniformScaling,n) = n*log(abs(x.λ))
factorize_(x) = factorize(x)
factorize_(x::UniformScaling) = x

function logpdf_mvt(v,μ,Σ,x)
  Σ_ = factorize_(Σ)
  d = length(x)
  C = (lgamma((v+d)/2)) - (lgamma(v/2)+(d/2)log(v*π)) - (0.5 * (logdet_(Σ_,d)))
  diff = abs.(μ.-x)
  # with poorly conditioned covariances, normdiff can be negative even when the
  # log determant is a valid, non-infinite number
  # NOTE: this somewhat akward notation is due to some missing method
  # definitions in the sparse matrix library (only right divide is defined)
  normdiff = max(0,(diff'*(Σ_'\diff))[1])
  result = C + -(v+d)/2*log(1+normdiff/v)
  @assert !isinf(result)
  result
end

mutable struct RidgeMultiNormalStats{T,C} <: Stats{T}
  μ::Vector{T}
  S::Vector{T}
  n::Float64
  x2_offset::Float64
  corr::C
end
Statistics.std(x::RidgeMultiNormalStats) = x.n > 0 ? x.S ./ x.n : Inf
Statistics.zero(x::RidgeMultiNormalStats{T},C::Coherence) where T =
  RidgeMultiNormalStats(zero(x.μ),zero(x.S),0.0,0.0,x.corr)

function findcorr(dims,dist,thresh)
  n = prod(dims)
  corr = spzeros(n,n)
  # TODO: consider using a ToeplitzMatrix
  # here
  for (i,ii) in enumerate(CartesianIndices(dims))
    for (j,jj) in enumerate(CartesianIndices(dims))
      if i <= j
        val = exp(-dist(ii.I,jj.I))
        if val > thresh
          corr[i,j] = val
        end
      end
    end
  end
  # @info "Sparsity of source-prior ($(size(corr,1)) × $(size(corr,2))) ridge: "*
  #   "$(nnz(corr)) ($(round(Int,100nnz(corr) / length(corr)))%)" maxlog=1

  Symmetric(corr)
end

function ridgenorm(prior::Coherence,n::Number,
                   x2_offset::Number=size(prior,2)*size(prior,3);
                   thresh=1e-3,scale=nothing,
                   freq=nothing)
  error("Outdated prior")
  ridgenorm(prior,n,rdist(scale=scale,freq=freq),x2_offset,thresh=thresh)
end

function ridgenorm(prior::Coherence,n::Number,dist::Function,
                   x2_offset::Number=size(prior,2)*size(prior,3);thresh=1e-3)
  error("Outdated prior")
  data = reshape(mean(prior,4),size(prior,1),:)
  μ = squeeze(mean(data,1),1)
  S = squeeze(sum(data.^2,1),1) .* n./size(data,1) .+ n.*thresh.^2
  corr = findcorr((size(prior,2),size(prior,3)),dist)
  RidgeMultiNormalStats{eltype(data)}(μ,S,n,x2_offset,corr)
end

# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
function update!(stats::RidgeMultiNormalStats{T},x::AbstractVector{T},
                 w=1.0) where T
  stats.n += w

  μ₀ = copy(stats.μ)
  stats.μ .+= (x .- stats.μ).*(w ./ stats.n)
  stats.S .+= (x .- stats.μ).*(x .- μ₀).*w

  stats
end

function downdate!(stats::RidgeMultiNormalStats{T},x::AbstractVector{T},
                   w=1.0) where T
  if w ≈ stats.n
    stats.μ .= 0
    stats.S .= 0
    stats.n = 0
  elseif w < stats.n
    μₙ = copy(stats.μ)
    stats.μ .= (stats.μ .- x.*(w ./ stats.n)) ./ (1 .- (w ./ stats.n))
    stats.S .-= (x .- stats.μ).*(x .- μₙ).*w

    stats.n -= w
  else
    error("Sum of weights would become negative. This happens when trying ",
          "to remove more samples then were added using `update!`.")
  end

  stats
end

function mult!(stats::RidgeMultiNormalStats,c)
  stats.n *= c

  stats
end

function Base.:(+)(a::RidgeMultiNormalStats{T},b::RidgeMultiNormalStats{T}) where T
  n = a.n + b.n
  RidgeMultiNormalStats(a.μ.*(a.n./n) .+ b.μ.*(b.n./n), a.S.+b.S,
                        n,a.x2_offset+b.x2_offset,a.corr)
end

function logpdf(stats::RidgeMultiNormalStats,x::AbstractVector)
  d = length(stats.μ)
  σ = sqrt.(stats.S ./ stats.n)
  nu = stats.n + stats.x2_offset - d + 1

  c = ((stats.n + 1)/(stats.n * nu))
  I,J,V = findnz(stats.corr.data)
  Λ = Symmetric(sparse(I,J,[v*σ[i]*σ[j] * c for (i,j,v) in zip(I,J,V)]))
  # Λ = ((Diagonal(σ) * stats.corr) * Diagonal(σ)) .* ((stats.n + 1)/(stats.n * nu))

  if nu < 1
    if length(x) != length(stats.μ)
      error("Dimension mistmatch between data ($(length(x))) and ",
            "the model ($(length(stats.μ))).")
    else
      error("Insufficient data for predictive distribution, ",
            "try collecting more data or increasing x2_offset")
    end
  end

  logpdf_mvt(nu,stats.μ,Λ,x)
end

struct ConstRidgePrior{T,C} <: Stats{T}
  S::T
  N::Float64
  x2_offset::Float64
  corr::C
end

function ridgenorm(prior::Number,N::Number,dims::Tuple,dist::Function,
                   x2_offset::Number=prod(dims);threshold)
  ConstRidgePrior(float(prior),Float64(N),Float64(x2_offset),
                  findcorr(dims,dist,threshold))
end

function ridgenorm(prior::Number,N::Number,dims::Tuple,
                   x2_offset::Number=prod(dims);scale,freq,threshold)
  dist = rdist(scale=scale,freq=freq)
  ridgenorm(float(prior),N,dims,dist,x2_offset,threshold=threshold)
end

function Base.zero(prior::ConstRidgePrior,d)
  RidgeMultiNormalStats(fill(zero(prior.S),d),fill(zero(prior.S),d),0.0,0.0,prior.corr)
end

function Base.:(+)(a::ConstRidgePrior{T},b::RidgeMultiNormalStats{T}) where T
  d = length(b.μ)
  n = b.n + a.N
  RidgeMultiNormalStats(b.μ.*(b.n/n),a.S.+b.S,n,a.x2_offset+b.x2_offset,
                        a.corr)
end

function logpdf(stats::ConstRidgePrior,x::AbstractVector)
  d = length(x)
  σ² = stats.S / stats.N
  Λ = stats.corr .* σ² .* ((stats.N + 1)/(stats.N * nu))
  nu = max(1,floor(Int,stats.N + stats.x2_offset - d + 1))

  logpdf_mvt(nu,0,Λ,x)
end
