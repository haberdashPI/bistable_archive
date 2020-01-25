using DSP

lowpass = digitalfilter(Lowpass(0.27),Butterworth(4))
function find_peaks(x)#,stage=:final)
  x = filtfilt(lowpass,x)
  dx = diff(x)
  ddx = diff(dx)
  cross0 = find(dx[1:end-1].*dx[2:end] .< 0) .+ 1

  # stage == :cross && return cross0

  peaks = cross0[ddx[cross0 .- 1] .< 0]

  # stage == :max && return peaks

  if !isempty(peaks)
    peaks[x[peaks] .> 0.75maximum(x[peaks])]
  else
    Int[]
  end
end

function guess_source_count(window)
  allpeaks = map(1:size(window,1)) do i
    window[i,find_peaks(window[i,:])]
  end

  maxi = map(allpeaks) do peaks
    isempty(peaks) ? typemin(eltype(peaks)) : maximum(peaks)
  end |> indmax
  length(allpeaks[maxi])
end

function slicerange(x,range,axis)
  dim = axisdim(x,axis)
  args = Array{Any}(ndims(x))
  args .= Colon()
  args[dim] = range
  getindex(x,args...)
end

function map_window(fn,x,window,delta;showprogress=false)
  indices = (1:floor(Int,delta / Δt(x)):ntimes(x))[1:end-1]
  width = floor(Int,window / Δt(x))
  vals = map(indices) do i
    fn(slicerange(x,clamp.(i:i+width,1,ntimes(x)),Axis{:time}))
  end
  AxisArray(vals,Axis{:time}(linspace(0s,duration(x),length(vals))))
end

function source_bumps(cs;window=1s,delta=0.25s)
  indices = 1:floor(Int,delta / Δt(cs)):ntimes(cs)
  y = similar(cs,Axis{:time}(times(cs)[indices]),AxisArrays.axes(cs)[2:end]...)

  for (k,i) in enumerate(indices)
    win = cs[atindex(0s .. window,i)]
    meanwin = meanabs(win,1)
    y[Axis{:time}(k)] = filtfilt(lowpass,copy(meanwin'))'
  end

  y
end

function source_peaks(cs;window=1s,delta=0.25s,stage=:final)
  indices = 1:floor(Int,delta / Δt(cs)):ntimes(cs)
  y = similar(cs,Axis{:time}(times(cs)[indices]),AxisArrays.axes(cs)[2:end]...)

  for (k,i) in enumerate(indices)
    win = cs[atindex(0s .. window,i)]
    meanwin = meanabs(win,1)
    peaks = mapslices(meanwin,2) do col
      ixs = find_peaks(col,stage)
      col .= 0
      col[ixs] = 1.0

      col
    end
    y[Axis{:time}(k)] = peaks
  end

  y
end

function source_count_by_peaks(cs;window=1s,delta=0.25s,buildup=1s)
  map_window(cs[buildup .. last(times(cs))],window,delta) do win
    guess_source_count(meanabs(win,1))
  end
end
