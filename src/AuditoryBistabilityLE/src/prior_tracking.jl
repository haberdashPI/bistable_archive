using Combinatorics
using DataStructures

struct PriorTracking{C,Sp,Fp} <: Tracking
  cohere::C
  tc::typeof(1.0s)
  source_prior::Sp
  freq_prior::Fp
  max_sources::Int
  min_norm::Float64
end
function checknorm(normalize,min_norm)
  if min_norm<Inf && normalize==false
    error("You must set normalize to true for `min_norm` value to be",
          " non-infinite")
  elseif min_norm == Inf && normalize==true
    error("You must specify a non-infinite value for `min_norm` if",
          " `normalize==true`.")
  end
end

Δt(x::PriorTracking) = Δt(x.cohere)
function Tracking(C,::Val{:prior};tc=1s,source_prior=nothing,
                  freq_prior=nothing, max_sources=4,
                  normalize=false,min_norm=Inf)
  checknorm(normalize,min_norm)
  PriorTracking(getmeta(C),tc,source_prior,freq_prior,max_sources,
                min_norm)
end

include("tracking_priors.jl")

# the number of ways observations A can be mapped to a set of source in B when
# it is possible that source b ∈ B is modeled as sum of elements in A
possible_groupings(n_sources,n_obs) =
  (Grouping(grouping,mapping)
   for grouping in partitions(1:n_obs)
   for sources in combinations(1:n_sources,length(grouping))
   for mapping in permutations(sources))

struct PermutedCoherence{C,P}
  cohere::C
  permuted::P
end

function prepare_coherence(C::Coherence,min_norm)
  @assert axisdim(C,Axis{:time}) == 1
  @assert axisdim(C,Axis{:scale}) == 2
  @assert axisdim(C,Axis{:freq}) == 3
  @assert axisdim(C,Axis{:component}) == 4

  # reorder for cache efficiency and then and normalize the data if necessary
  C_ = permutedims(C,[2,3,4,1])
  if min_norm < Inf
    for I in CartesianIndices((Base.axes(C_,3),Base.axes(C_,4)))
      C_[:,:,I] ./= max(min_norm,norm(vec(view(C_,:,:,I))))
    end
  end
  PermutedCoherence(C,C_)
end

function track(C::Coherence,params,args...) 
  track(prepare_coherence(C,params.min_norm),params,args...)
end
function track(perm::PermutedCoherence,params::PriorTracking,progressbar=true,
               progress = track_progress(progressbar,ntimes(C),"prior"))
  C = perm.cohere
  C_ = perm.permuted
  timeax = 4

  track = TrackedSources(size(C_,1)*size(C_,2),params)

  C_out = zeros(eltype(C_),size(C_,1),size(C_,2)...,
                params.max_sources,size(C_,4))
  source_out = copy(C_out)
  sourceS_out = copy(C_out)

  window = ceil(Int,params.tc / Δt(C))
  lp_out = AxisArray(fill(0.0,ntimes(C)),AxisArrays.axes(C,1))
  buf = Array{eltype(C_)}(undef,size(C_,1),size(C_,2))
  function sumcomponents(C_t,kk)
    y = sum!(buf,view(C_t,:,:,kk,:))
    # buf ./= size(C_t,timeax)
    vec(buf)
  end

  veccomponent(C_t,k) = vec(view(C_t,:,:,k))

  groupings = CircularDeque{Grouping}(window+1)
  for t in eachindex(times(C))
    # find the MAP grouping of sources
    MAPgrouping, lp_out[t] = 
      maximumby(g -> logpdf(track,view(C_,:,:,:,t:t),sumcomponents,g),
                possible_groupings(params.max_sources,ncomponents(C)))

    # arrange the sources
    for (kk,i) in iterable(MAPgrouping)
      for k in kk
        C_out[:,:,i,t] .+= C_[:,:,k,t]  
      end
    end

    # update the source models
    update!(track,view(C_out,:,:,:,t),veccomponent,MAPgrouping)
    # mult!(track,decay)
    push!(groupings,MAPgrouping)
    if t > window
      downdate!(track,view(C_out,:,:,:,t-window),veccomponent,popfirst!(groupings))
    end

    next!(progress)
  end


  tracking =
    MetaArray(params,AxisArray(C_out,AxisArrays.axes(C,2),
                                     AxisArrays.axes(C,3),
                                     Axis{:component}(1:params.max_sources),
                                     AxisArrays.axes(C,1)))
  (tracking,lp_out)
end

function sort_components(x)
  order = sortperm(component_means(x),rev=true) 
  x .= x[Axis{:component}(order)]
  x
end
