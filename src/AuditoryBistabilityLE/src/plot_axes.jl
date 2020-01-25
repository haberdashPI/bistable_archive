using PlotAxes

function PlotAxes.asplotable(x::AxisArray{<:SourceTracking},args...;component=1,kwds...)
  vals = cat((xi[:,:,component,:] for xi in x)...,dims=4)
  axs = AxisArrays.axes(x[1])
  start,stop = extrema(log.(ustrip.(freqs(x[1]))))
  logfreq = Axis{:logfreq}(range(start,stop=stop,length=length(freqs(x[1]))))
  y = AxisArray(vals,axs[1],logfreq,axs[4],AxisArrays.axes(x)...)
  asplotable(y,:time,:logfreq,:scale,:params,args...;kwds...)
end
