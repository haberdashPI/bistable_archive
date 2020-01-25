using MacroTools
using Match
using Parameters
using Parameters
using AxisArrays
using Unitful

export cohere, component, mask, ncomponents, components, component_means,
  windowing, map_windowing

abstract type CoherenceMethod end
struct CParams{M,P}
  cort::P
  ncomponents::Int
  skipframes::Int
  window::typeof(1.0s)
  delta::typeof(1.0s)
  method::M
end

const Coherence = MetaArray{<:AxisArray,<:CParams}
const CParamLike = Union{CParams,Coherence}

CorticalSpectralTemporalResponses.Δt(as::CParamLike) = as.delta

function Base.show(io::IO,::MIME"text/plain",x::Coherence)
  if !get(io, :compact, false)
    println(io,"Coherence Components")
    describe_axes(io,x)
  else
    println(io,string(duration(x))," Coherence Components")
  end
end

ncomponents(x::CParams) = x.ncomponents
ncomponents(x::AbstractArray) = size(x,axisdim(x,Axis{:component}))
components(x::CParams) = 1:ncomponents(x)
components(x::AbstractArray) = axisvalues(AxisArrays.axes(x,Axis{:component}))[1]
component(x::AbstractArray,n) = x[Axis{:component}(n)]

CorticalSpectralTemporalResponses.hastimes(x::Coherence) = HasTimes()

function component_means(C)
  mdims = filter(x -> x != axisdim(C,Axis{:component}),1:ndims(C))
  vec(mean(C,dims=mdims))
end

CorticalSpectralTemporalResponses.frame_length(params::CParams,x) =
  max(1,floor(Int,params.delta / Δt(x)))

function CParams(x;ncomponents=1,window=1000ms,
                 delta=10ms,method=:nmf,skipframes=0,method_kwds...)
  method = CoherenceMethod(Val{method},method_kwds)

  CParams(getmeta(x), ncomponents, skipframes,
          convert(typeof(1.0s),asseconds(window)), convert(typeof(1.0s),asseconds(delta)), method)
end

windowing(x,params::CParams) =
  windowing(x,length=params.window,step=params.delta)
windowlen(params::CParams,x) = round(Int,params.window/Δt(x))

function nunits(params::CParams,x)
  mapreduce(*,AxisArrays.axes(x)) do ax
    isa(ax,Axis{:time}) || isa(ax,Axis{:rate}) ? 1 : length(ax)
  end
end

cohere(x::MetaUnion{AxisArray};progressbar=true,params...) =
  cohere(x,CParams(x;params...),progressbar)

function cohere_progress(progressbar,x,params)
  if progressbar
    windows = windowing(x,params)
    Progress(length(windows),desc="Temporal Coherence Analysis: ")
  end
end

function cohere(x::MetaUnion{AxisArray},params::CParams,
                progressbar=true,
                progress = cohere_progress(progressbar,x,params))
  @assert axisdim(x,Axis{:time}) == 1
  @assert axisdim(x,Axis{:rate}) == 2
  @assert axisdim(x,Axis{:scale}) == 3
  @assert axisdim(x,Axis{:freq}) == 4
  @assert ndims(x) == 4

  # if we already have components, just wrap up the values with
  # parameters (since we've already computed components)
  if :component in axisnames(x)
    return MetaArray(params,x)
  end

  windows = windowing(x,params)

  K = ncomponents(params)
  C_data = zeros(eltype(params.method,x),length(windows),size(x)[3:end]...,K)
  C = AxisArray(C_data,
                AxisArrays.axes(windows,Axis{:time}),
                AxisArrays.axes(x)[3:end]...,
                Axis{:component}(1:K))

  with_method(params.method,K) do extract
    for (i,w_inds) in enumerate(windows)
      skipped = w_inds[1:1+params.skipframes:end]
      # axis array can't handle skipped indices, so we assume
      # the right dimensionality
      components = extract(x[skipped,:,:,:])
      C[i,Base.axes(components)...] = components

      next!(progress)
    end
  end

  MetaArray(params,C)
end
