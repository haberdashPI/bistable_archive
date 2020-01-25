
using AxisArrays

windowing(x,dim=timedim(x);kwds...) = windowing(hastimes(x),x,dim;kwds...)
map_windowing(fn,x,dim=timedim(x);kwds...) =
  map_windowing(fn,hastimes(x),x,dim;kwds...)

function map_windowing(fn,::HasTimes,x,dim;step=nothing,flatten=false,kwds...)
  windows = windowing(x,dim;step=step,kwds...)
  xs = map(windows) do ixs
    fn(x[Axis{:time}(ixs)])
  end

  if flatten
    result = cat(xs...,dims = ndims(xs[1])+1)
    axisnames = Symbol.("ax".*string.(1:ndims(result)-1))
    AxisArray(result,
              (Axis{ax}(1:n) for (ax,n) in zip(axisnames,
                                               size(result)[1:end-1]))...,
              AxisArrays.axes(windows,Axis{:time}))
  else
    AxisArray(xs,AxisArrays.axes(windows,Axis{:time}))
  end
end
map_windowing(fn,::HasNoTimes,x,dim;kwds...) =
  map(ixs -> fn(x[Axis{:time}(ixs)]),windowing(HasNoTimes(),x,dim;kwds...))

function windowing(::HasNoTimes,x::AbstractArray,dim;
                   length=nothing,step=nothing,minlength=length)
  (max(1,t-length+1):t for t in Base.axes(x,dim)[min(minlength,end):step:end])
end

function windowing(::HasTimes,data::AbstractArray,dim;
                   length=nothing,step=nothing,minlength=length)
  helper(x::Number) = max(1,floor(Int,x*s / Δt(data)))
  helper(x::Quantity) = max(1,floor(Int,x / Δt(data)))
  length_,step_,minlength_ = helper.((length,step,minlength))

  win = windowing(HasNoTimes(),data,dim,
                  length=length_,step=step_,minlength=minlength_)
  AxisArray(collect(win),
            Axis{:time}((min(minlength_,size(data,dim)):step_:size(data,dim)).*
                        Δt(data) .- (length_/2*Δt(data))))
end
