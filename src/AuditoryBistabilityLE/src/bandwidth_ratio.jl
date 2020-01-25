export bandwidth_ratio

function fixlength(x,n;fill=0)
  if length(x) > n
    view(x,1:n)
  elseif length(x) < n
    append!(convert(Array{Union{eltype(x),Missing}},x),
            Base.fill(fill,n-length(x)))
  else
    x
  end
end

function dominant_channels(sp;level_threshold=0.9,threshold=0.25,
                           window=500ms,delta=250ms,thresh_ratio=1.0,
                            max_channels=10)
  x = map_windowing(sp,length=window,step=delta,flatten=true) do window
    thresh_level = thresh_ratio*quantile(vec(sp),1 - threshold/size(sp,2))
    # @show thresh_level
    levels = dropdims(mapslices(x -> quantile(x,level_threshold),
                                Array(window),dims=1), dims=1)
    # @show sort(levels,rev=true)[1:min(end,20)]
    # thresh_level = thresh_ratio*
    #   quantile(levels,1 - threshold/length(levels))

    x = fixlength(findall(levels .> thresh_level),max_channels,fill=missing)
    # @show x
    x
  end 
end

function bandwidth_ratio(spmask::AbstractMatrix, sp::AbstractMatrix,
                         settings)
  settings = read_settings(settings)
  startHz, stopHz = asHz.(settings.rates.freqbound.freq_limits)
  bandwidth_ratio(spmask, sp[:,startHz .. stopHz]; settings.bandwidth_ratio...)
end

function bandwidth_ratio(spmask, sp; threshold=1.5,
                         window=500ms,
                         full_band_ratio=2,
                         level_threshold=0.9,
                         thresh_ratio=1.0,
                         delta=250ms)
  window = asseconds(window)
  @assert frequencies(spmask) == frequencies(sp) "Frequency axes are not equal."
  @assert times(spmask) == times(sp) "Time axes are not equal."

  fullchan = dominant_channels(sp,window=window*full_band_ratio,
                               delta=delta,threshold=threshold,
                               level_threshold=level_threshold,
                               thresh_ratio=thresh_ratio)

  maskchan = dominant_channels(spmask,window=window,delta=delta,
                               threshold=threshold,
                               level_threshold=level_threshold,
                               thresh_ratio=thresh_ratio)

  dim = axisdim(fullchan,Axis{:time})
  others = setdiff(1:ndims(fullchan),[dim])
  fullband = mapslices(Array(fullchan),dims=others) do x
    maximum(skipmissing(x)) - minimum(skipmissing(x))
  end |> x -> AxisArray(vec(x),Axis{:time}(times(fullchan)))

  maskband = AxisArray(Array{Float64}(undef,size(maskchan,dim)),
                       Axis{:time}(times(maskchan)))
  ratio = AxisArray(Array{Float64}(undef,size(maskchan,dim)),
                    Axis{:time}(times(maskchan)))
  winwidth = window*full_band_ratio
  for (i,t) in enumerate(times(maskchan))
    val = (-window*full_band_ratio*0.51 .. window*full_band_ratio*0.51) + t
    maskwin = view(maskchan,Axis{:time}(t))
    fullwin = view(fullchan,Axis{:time}(val))
    maskcount = skipmissing(intersect(maskwin,fullwin))
    if isempty(maskcount)
      maskband[i] = 0
    else
      maskband[i] = maximum(maskcount) - minimum(maskcount)
    end
    ratio[i] = maskband[i] / maximum(view(fullband,Axis{:time}(val)))
  end

  ratio, fullband, maskband
end

