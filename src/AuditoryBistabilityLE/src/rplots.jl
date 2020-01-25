using .RCall
# import CorticalSpectralTemporalResponses: rplot, raster_plot
export rplot, scale_plot, cohere_plot

R"""
library(ggplot2)
"""

function titlefn(λ,λ_digits,components)
  @assert !any(isnan,λ)
  digits = λ_digits != :automatic ?  λ_digits :
    min(floor(Int,-log10(clamp(λ[minimum(components)],1e-8,1e8)))+2,
        floor(Int,-log10(clamp(λ[maximum(components)],1e-8,1e8)))+1)

  function title(n)
    if n <= length(λ)
      # nstr = @sprintf("%02d",n)
      "lambda[$(n)] == $(round(λ[n],digits=digits))"
      # "Lmb_$nstr = $(round(λ[n],3))"
    else
      "sigma^2 == $(round(sum(var(C)),digits=digits))"
      # "Variance = $(sum(var(C)))"
    end
  end
end

function nearin(xin,xs)
  _,inds = findmin(abs.(vec(xin) .- vec(xs)'),2)
  cols = map(ii -> ii[2],CartesianIndices((length(xin),length(xs))))
  if length(inds) > 0
    cols[inds[:,1]]
  else
    Array{Int}(0)
  end
end

function findnear(x,nearby)
  indices = sort(unique(nearin(filter(!isnan,nearby),
                                filter(!isnan,x))))
  if any(isnan.(nearby))
    [NaN; x[indices]], [find(isnan,x); indices]
  else
    x[indices], indices
  end
end

function cohere_plot(slice::AxisArray{T,3} where T; kwds...)
  @assert axisdim(slice,Axis{:scale}) == 1
  @assert axisdim(slice,Axis{:freq}) == 2
  @assert axisdim(slice,Axis{:component}) == 3

  ii = CartesianIndices(size(slice))
  at(i) = map(ii -> ii[i],ii)

  df = DataFrame(value = vec(slice),
                 scale = vec(round.(ustrip.(uconvert.(cycoct,
                                                      scales(slice)[at(1)])),
                                    digits=2)),
                 freq_bin = vec(at(2)),
                 component = vec(at(3)))

  fbreaks,findices = freq_ticks(slice)
  p = raster_plot(df;value=:value,x=:scale,y=:freq_bin,kwds...)

R"""
  $p + facet_grid(component~.) + #ordered_scales(scale)) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    ylab('Frequency (Hz)') + xlab('Time (s)')
"""
end

# function rplot(C::Coherence{<:Any,<:Any,<:Any,4} where {M,T};λ_digits=:automatic,
#                kwds...)
#   C = sort_components(C)
#   @assert axisdim(C,Axis{:time}) == 1
#   @assert axisdim(C,Axis{:scale}) == 2
#   @assert axisdim(C,Axis{:freq}) == 3
#   @assert axisdim(C,Axis{:component}) == 4

#   ii = CartesianIndices(size(C))
#   at(i) = map(ii -> ii[i],ii)

#   rowtitle = titlefn(component_means(C),λ_digits,Base.axes(C,4))

#   df = DataFrame(value = vec(C),
#                  time = vec(ustrip.(uconvert.(s,times(C)[at(1)]))),
#                  scale = vec(round.(ustrip.(uconvert.(cycoct,scales(C)[at(2)])),
#                                     digits=2)),
#                  freq_bin = vec(at(3)),
#                  component = vec(at(4)),
#                  component_title = rowtitle.(vec(at(4))))

#   fbreaks,findices = freq_ticks(C)
#   p = raster_plot(df;value=:value,x=:time,y=:freq_bin,kwds...)

# R"""
#   scalevals = $(ustrip.(uconvert.(cycoct,scales(C))))
#   scalestr = function(x){sprintf("'Scale: %3.2f cyc/oct'",x)}
#   ordered_scales = function(x){
#     factor(scalestr(x),levels=scalestr(scalevals))
#   }

#   $p + facet_grid(ordered_scales(scale)~component_title,labeller=label_parsed) +
#   scale_y_continuous(breaks=$findices,labels=$fbreaks) +
#   ylab('Frequency (Hz)') + xlab('Time (s)')
# """
# end

function rplot(S::SourceTracking;λ_digits=:automatic,
               kwds...)
  S = sort_components(S)
  @assert axisdim(S,Axis{:scale}) == 1
  @assert axisdim(S,Axis{:freq}) == 2
  @assert axisdim(S,Axis{:component}) == 3
  @assert axisdim(S,Axis{:time}) == 4

  ii = CartesianIndices(size(S))
  at(i) = map(ii -> ii[i],ii)

  rowtitle = titlefn(component_means(S),λ_digits,Base.axes(S,3))

  df = DataFrame(value = vec(S),
                 time = vec(ustrip.(uconvert.(s,times(S)[at(4)]))),
                 scale = vec(round.(ustrip.(uconvert.(cycoct,scales(S)[at(1)])),
                                    digits=2)),
                 freq_bin = vec(at(2)),
                 component = vec(at(3)),
                 component_title = rowtitle.(vec(at(3))))

  fbreaks,findices = freq_ticks(S)
  p = raster_plot(df;value=:value,x=:time,y=:freq_bin,kwds...)

R"""
  scalevals = $(ustrip.(uconvert.(cycoct,scales(S))))
  scalestr = function(x){sprintf("'Scale: %3.2f cyc/oct'",x)}
  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr(scalevals))
  }

  $p + facet_grid(ordered_scales(scale)~component_title,labeller=label_parsed) +
  scale_y_continuous(breaks=$findices,labels=$fbreaks) +
  ylab('Frequency (Hz)') + xlab('Time (s)')
"""
end


function rplot(tracks::AxisArray{<:SourceTracking},component=1;kwds...)
  data = AxisArray(cat((sort_components(track)[:,:,component,:]
                        for track in tracks)...,dims=4),
                   AxisArrays.axes(tracks[1])[[1,2,4]]...,
                   Axis{:track}(1:length(tracks)))

  ii = CartesianIndices(size(data))
  at(i) = vec(map(ii -> ii[i],ii))
  track_vals = axisvalues(tracks)[1]

  df = DataFrame(value = vec(data),
                 time = vec(ustrip.(uconvert.(s,times(data)[at(3)]))),
                 scale = vec(round.(ustrip.(uconvert.(cycoct,
                                                      scales(data)[at(1)])),
                                    digits=2)),
                 freq_bin = vec(at(2)),
                 tau = vec([ustrip(τ) for (τ,_,__) in track_vals[at(4)]]),
                 sigma = vec([ustrip(σ) for (_,σ,__) in track_vals[at(4)]]),
                 N = vec([N for (_,__,N) in track_vals[at(4)]]))

  # @show(df[1:5,:])

  fbreaks,findices = freq_ticks(tracks[1])
  p = raster_plot(df;value=:value,x=:time,y=:freq_bin,kwds...)

R"""
  scalevals = $(ustrip.(uconvert.(cycoct,scales(data))))
  scalestr = function(x){sprintf("'Scale: %3.2f cyc/oct'",x)}
  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr(scalevals))
  }

  $p + facet_grid(ordered_scales(scale)~tau+sigma+N,labeller=label_both) +
  scale_y_continuous(breaks=$findices,labels=$fbreaks) +
  ylab('Frequency (Hz)') + xlab('Time (s)')
"""
end

function rplot(tracks::AxisArray{<:Coherence},component=1;kwds...)
  data = AxisArray(cat(4,(track[:,:,:,component] for track in tracks)...),
                   AxisArrays.axes(tracks[1])[1:3]...,Axis{:track}(1:length(tracks)))

  ii = CartesianIndices(size(data))
  at(i) = vec(map(ii -> ii[i],ii))
  track_vals = axisvalues(tracks)[1]
  # @show track_vals
  # @show unique(at(4))

  df = DataFrame(value = vec(data),
                 time = vec(ustrip.(uconvert.(s,times(data)[at(1)]))),
                 scale = vec(round.(ustrip.(uconvert.(cycoct,
                                                      scales(data)[at(2)])),
                                    digits=2)),
                 freq_bin = vec(at(3)),
                 tau = vec([ustrip(τ) for (τ,_) in track_vals[at(4)]]),
                 sigma = vec([ustrip(σ) for (_,σ) in track_vals[at(4)]]))
  # @show(df[1:5,:])

  fbreaks,findices = freq_ticks(tracks[1])
  p = raster_plot(df;value=:value,x=:time,y=:freq_bin,kwds...)

R"""
  scalevals = $(ustrip.(uconvert.(cycoct,scales(data))))
  scalestr = function(x){sprintf("'Scale: %3.2f cyc/oct'",x)}
  ordered_scales = function(x){
    factor(scalestr(x),levels=scalestr(scalevals))
  }

  $p + facet_grid(ordered_scales(scale)~tau+sigma,labeller=label_parsed) +
  scale_y_continuous(breaks=$findices,labels=$fbreaks) +
  ylab('Frequency (Hz)') + xlab('Time (s)')
"""
end

#=
function rplot(cohere::CoherenceModel,C::Factors;
               components=1:ncomponents(C),λ_digits=:automatic)
  u = factors(C)[:,components]
  u = reshape(u,length(scales(cohere)),:,size(u,2))

  ii = CartesianIndices(size(u))
  at(i) = vec(map(ii -> ii[i],ii))

  title = titlefn(strengths(C),λ_digits,components)

  df = DataFrame(response = vec(u),
                 scale_index = at(1),
                 freq_bin = at(2),
                 component = at(3),
                 component_title = title.(at(3)))

  sindices = 1:2:length(scales(cohere))
  sbreaks = round.(scales(cohere)[sindices],2)
  fbreaks,findices = freq_ticks(cohere)

  p = raster_plot(df,value=:response,x=:scale_index,y=:freq_bin)

R"""

  library(ggplot2)

  $p + facet_wrap(~component_title,labeller=label_parsed) +
    scale_y_continuous(breaks=$findices,labels=$fbreaks) +
    scale_x_continuous(breaks=$sindices,labels=$sbreaks) +
    ylab('Frequency (Hz)') + xlab('Scale (cycles/octave)')
"""

end

function rplot(tempc::CoherenceModel,C::AbstractArray{<:FactorSeries},
               label=:index => 1:length(C))
  x = [fusion_signal(tempc,Ci) for Ci in C]
  df = DataFrame(resp = vcat((real.(xi) for xi in x)...),
                 var = vcat((fill(label[2][i],length(x[i]))
                             for i in eachindex(x))...),
                             time = vcat((ustrip.(uconvert.(s,eachindex(xi) * Δt(tempc)))
                              for xi in x)...))

R"""
  library(ggplot2)

  ggplot($df,aes(x=time,y=resp,group=factor(var),color=factor(var))) +
    geom_line() +
    scale_color_brewer(palette='Set1',name=$(string(label[1]))) +
    coord_cartesian(ylim=c(1.1,0)) + ylab('lambda / var(x)') +
    xlab('time (s)')
"""
end


function strength_plot(cohere::CoherenceModel,C::EigenSeries;n=ncomponents(C))
  λ = eigvals(C)
  λ = λ[:,sortperm(abs.(λ[max(1,end-10),:]),rev=true)]
  ii = CartesianIndices(size(λ))
  at(i) = map(ii -> ii[i],ii)

  prop_complex = mean(imag.(λ))
  if prop_complex > 0
    warn("$(round(100prop_complex,1))% of eigenvalues were "*
         "complex (using absolute value).")
  end

  df = DataFrame(value = vec(abs.(λ) ./ sum.(var.(C))),
  time = vec(ustrip(uconvert.(s,at(1) * Δt(cohere)))),
                 component = vec(at(2)))

  df = df[df[:component] .<= n,:]

R"""
  ggplot($df,aes(x=time,y=value,color=factor(component),group=component)) +
    geom_line() + scale_color_brewer(palette='Set1',name='Component') +
    xlab('Time (s)') + ylab('Value')
"""
end

function rplot(cohere::CoherenceModel,λ::Vector)
  @assert all(imag.(λ) .== 0) "Can't plot complex eigenvalues"
  λ = real.(λ)
  df = DataFrame(value = sort(λ,rev=true),index = collect(eachindex(λ)))

R"""
  ggplot($df,aes(x=index,y=value)) + geom_bar(stat='identity')
"""
end
=#
