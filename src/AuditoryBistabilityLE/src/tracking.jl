export track, topN
using Combinatorics
using DataStructures

abstract type Tracking end
@with_kw struct SimpleTracking <: Tracking
  tc::typeof(1.0s) = 1s
end
Tracking(::Val{:simple};params...) = SimpleTracking(;params...)

# TODO: this makes me realize I can make a more general metadata array that
# isn't an axis array at all
const WithAxes{Ax} = AxisArray{<:Any,<:Any,<:Any,Ax}
const SourceTracking = 
  MetaArray{<:WithAxes{<:Tuple{Axis{:scale},Axis{:freq},
                               Axis{:component},Axis{:time}}}}
component_means(S::SourceTracking) = vec(mean(S,dims=[1,2,4]))
function component_means(S::SubArray{<:Any,<:Any,<:SourceTracking}) 
  vec(mean(S,dims=[1,2,4]))
end

function track(C::Coherence;method=:simple,progressbar=true,params...)
  method = Tracking(C,Val{method}();params...)
  track(C,method,progressbar)
end

function topN(by,N,xs)
  state = start(xs)
  if done(xs,state)
    error("Empty iterable.")
  end
  x, state = next(xs,state)
  by_x = by(x)
  queue = PriorityQueue{typeof(by_x),typeof(x)}()
  enqueue!(queue,x,by_x)

  while !done(xs,state)
    x, state = next(xs,state)
    by_x = by(x)
    enqueue!(queue,x,by_x)
    if length(queue) > N
      dequeue!(queue)
    end
  end

  result = Array{eltype(xs)}(length(queue))
  i = 0
  while !isempty(queue)
    result[i += 1] = dequeue!(queue)
  end
  result
end

function maximumby(by,xs)
  itr = iterate(xs)
  if itr == nothing
    error("Empty iterator.")
  end

  x,state = itr
  result = x
  maxval = by(x)
  itr = iterate(xs,state)
  while itr !== nothing
    x, state = itr
    val = by(x)
    if isless(maxval,val)
      result = x
      maxval = val
    end
    itr = iterate(xs,state)
  end
  result, maxval
end

function bestordering(x,y)
  @assert length(x) == length(y)
  @assert length(x) <= 6 "No more than 6 components supported"*
    " (due to combinatorial explosion of possible solutions.)"
  # minimize cosine difference
  diff(x,y) = norm(vec(x) .- vec(y),1)

  maximumby(permutations(1:length(y))) do order
    -sum(diff(xi,yi) for (xi,yi) in zip(x,y[order]))
  end
end

function track_progress(progressbar,n,name)
  if progressbar
    Progress(n,desc="Source Tracking ($name): ")
  end
end
