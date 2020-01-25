using Gadfly
using CSV
using DependentBootstrap
using DataFramesMeta
using ShiftedArrays
using PlotAxes
PlotAxes.set_backend!(:gadfly)

function packing(x;maxpad=true,pad=0.5)
    vals = sort!(unique(x))
    pos = [0; cumsum([1.0+pad; fill(1,length(vals)-3); maxpad ? 1.0+pad : 1])]
    vals, pos
end

function packaxes(x;kwds...)
    p = Dict(zip(packing(x;kwds...)...))
    [p[xi] for xi in x]
end

function packaxes_invfn(x;kwds...)
   p,v = packing(x;kwds...)
   p = Dict(zip(v,p))
   xi -> p[xi]
end

function rename_levels_for(df,vals;suffixes=[:c_a,:c_m])
  for suf in suffixes
    if Symbol("f_"*string(suf)) in names(df)
      df[suf] = zero(df[Symbol("f_"*string(suf))])
    elseif Symbol("s_"*string(suf)) in names(df)
      df[suf] = zero(df[Symbol("s_"*string(suf))])
    elseif Symbol("t_"*string(suf)) in names(df)
      df[suf] = zero(df[Symbol("t_"*string(suf))])
    else
      error("No level parameter with suffix $suf")
    end
  end
  df[:level] = "unknown"
  for (prefix,level) in [("f_","Peripheral"),("s_","Cortical"),("t_","Object")]
    rows = df[Symbol(prefix*"c_σ")] .> 0
    df[rows,:level] = level
    for suffix in suffixes
      if Symbol(prefix*string(suffix)) in names(df)
        df[rows,suffix] = df[rows,Symbol(prefix*string(suffix))]
      else
        df[rows,suffix] = NaN
      end
    end
  end
  df[[suffixes;:level;vals]]
end

function DataFramesMeta.linq(::DataFramesMeta.SymbolParameter{:rename_levels},
                             df, vals)
  :(rename_levels($df,$vals))
end

function colorscale(smap;reverse=false,minvalue,maxvalue,colorstart=minvalue,
                    colorstop=maxvalue,
                    midvalue=(colorstop-colorstart)/2+colorstart)
    scale(x,min,max) = (x - min)/(max - min)
    scalei(x,min,max) = x*(max - min) + min
    fi(x) = clamp(scalei(x,minvalue,maxvalue),minvalue,maxvalue)
    c(x) = scale(x,colorstart,colorstop)
    colors = !reverse ? colormap(smap,mid=c(midvalue)) :
      reverse!(colormap(smap,mid=1-c(midvalue)))

    function(x)
        colors[1+clamp(floor(Int,99*c(fi(x))),0,99)]
    end
end

function select_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                     stop_time=20s,kwds...)
  selections = select_params(params;Δf=Δf,kwds...)
  if length(selections) > 1
    error("Ambiguous parameter selection, mutliple matches.")
  end
  selection = selections[1]
  masks = []
  for_results_in(joinpath(datadir,"data")) do entry
    if entry["pindex"] == selection[1]
      push!(masks,entry["mask"])
    end
  end
  audiospect(masks[simulation],settings)[start_time .. stop_time]
end

function unzip(xs)
  fns = map(enumerate(first(xs))) do (i,el)
    ys -> map(y -> y[i],ys)
  end |> Tuple
  map(fns) do fn
    fn(xs)
  end
end

function qqdata(x,y)
  nx = length(x)
  ny = length(y)
  order = sortperm([x; y])
  qx,qy = mapreduce((x,a) -> push!(x,x[end].+a),order,init=[(0.0,0.0)]) do i
    i <= nx ? (1/nx,0.0) : (0.0,1/ny)
  end |> unzip

  qx, qy, order
end

function cumreduce(fn,xs)
  reduce(xs[2:end];init=[first(xs)]) do result,x
    push!(result,fn(result[end],x))
  end
end

function plot_stream_data(df,params,selections::Vector;exclude_human=false)
  sims = map(enumerate(selections)) do (i,sel)
    df_,params_ = select_data(df,params;sel...)
    sim = by(stream_summary(df_,params_),:st) do g
      DataFrame([findci(g.streaming)])
    end
    sim[!,:experiment] .= "simulation"
    sim[!,:sid] .= "sel"*string(i)

    sim
  end |> x -> vcat(x...)

  sim = by(sims,:st) do g
    DataFrame([findci(g.mean)])
  end
  sim[!,:experiment] .= "simulation"

  if !exclude_human
    hdata = human_data()
    human = by(hdata.stream,:st) do g
      DataFrame([findci(collect(skipmissing(g.streaming)))])
    end
    human[!,:experiment] .= "human"

    vcat(sim,human)
  else
    sim
  end
end

function plot_stream_data(df,params;exclude_human=false,kwds...)
  df,params = select_data(df,params;kwds...)
  hdata = human_data()

  sim = by(stream_summary(df,params),:st) do g
    DataFrame([findci(g.streaming)])
  end
  sim.experiment .= "simulation"

  if !exclude_human
    human = by(hdata.stream,:st) do g
      DataFrame([findci(g.streaming)])
    end
    human.experiment .= "human"

    vcat(sim,human)
  else
    sim
  end
end

function plot_stream(df,others...;kwds...)
  plot_stream_data(df,others...;kwds...) |> plot_stream
end

function plot_stream(stream)
  stream = @transform(stream,
                      pos = @.(log2(:st) -
                               ifelse(:experiment == "human",0.05,-0.05)));

  plot(stream,x=:pos,y=:mean,ymin=:lowerc,ymax=:upperc,
     color=:experiment,shape=:experiment,
     Guide.shapekey(pos=[log2(3),1.4]),
     Guide.xticks(ticks=log2.([3,6,12])),
     Scale.color_discrete_manual("lightgray","darkgray"),
     Scale.x_continuous(labels=x -> string(floor(Int,2.0^x))),
     Scale.y_continuous(labels=x -> string(floor(Int,100*x))),
     Guide.xlabel("Δf (semitones)",orientation=:horizontal),
     Guide.ylabel("% streaming",orientation=:vertical),
     Coord.cartesian(xmin=log2(3)-0.25,xmax=log2(12)+0.25),
     Geom.point,Geom.line,Geom.errorbar,
     Theme(discrete_highlight_color=x->"black"))
end

# if a line doesn't move stepwise (always horizontal or vertical lines) add
# points in-between so that it does move stepwise
function addsteps(xs,ys)
  xr = eltype(xs)[]
  yr = eltype(ys)[]
  oldx = 0
  oldy = 0
  for (x,y) in zip(xs,ys)
    if oldx != x
      push!(xr,x)
      push!(yr,oldy)
    end
    push!(xr,x)
    push!(yr,y)

    oldx = x
    oldy = y
  end
  xr, yr
end

function cleanline(xs,ys)
  indices = filter(2:length(xs)-1) do i
    !(xs[i-1] == xs[i] == xs[i+1]) ||
    !(ys[i-1] == ys[i] == ys[i+1])
  end
  xs[[1;indices;end]], ys[[1;indices;end]]
end

function plot_lengths_data(df,params,selections::Vector;norm=identity,
  minlength=0.5,kwds...)

  lengths = map(enumerate(selections)) do (i,sel)
    df_,params_ = select_data(df,params;sel...)
    lens = length_summary(df_,params_)
    DataFrame(lengths = lens,sid = i)
  end |> x -> vcat(x...)

  sum_prop = 0.0
  normed = by(lengths,:sid) do df
    sum_prop += mean(df.lengths .< minlength)
    lens = filter_minlength(minlength,string("sim-",first(df.sid)),
      df.lengths)
    DataFrame(sid = first(df.sid), lengths = norm(lens))
  end
  slens = normed.lengths

  mean_percent = round(100sum_prop / length(unique(lengths.sid)),digits=2)
  @info "The average simulation response less than minlength=$minlength is:"*
    " %$mean_percent."

  hdata = human_data()
  humanlens = hdata.lengths
  sum_prop = 0.0
  normed = by(humanlens,:sid) do df
    sum_prop += mean(df.lengths .< minlength)
    lens = filter_minlength(minlength,string("human-",first(df.sid)),
      df.lengths)
    DataFrame(sid = first(df.sid), lengths = norm(lens))
  end
  hlens = normed.lengths

  mean_percent = round(100sum_prop / length(unique(humanlens.sid)),digits=2)
  @info "The average human response less than minlength=$minlength is:"*
    " %$mean_percent."

  (human=hlens,simulation=slens)
end

function plot_lengths_qq(df,others...;kwds...)
  plot_lengths(plot_lengths_data(df,others...);kwds...)
end

function plot_lengths_hist((hlens,slens)::NamedTuple;xmax=20,binstep=1,
                           density=false)
  if density
    df = DataFrame(length = [hlens;slens],
                   experiment = [fill("human",length(hlens));
                                 fill("simulation",length(slens))])
    plot(df,x=:length,color=:experiment,
         Stat.density, Geom.polygon(fill=true,preserve_order=true),
         Guide.colorkey(pos=[0.65*Gadfly.w,-0.3*Gadfly.h]),
         Theme(lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.4)),
         Coord.cartesian(xmax=xmax))
  else
    if any(hlens .> xmax)
      @warn("Ignoring ~$(100round(mean(hlens .> xmax)))% of human data")
    elseif any(slens .> xmax)
      @warn("Ignoring ~$(100round(mean(slens .> xmax)))% of simulation data")
    end
    hheight = normalize(fit(Histogram,hlens,0:binstep:xmax))
    sheight = normalize(fit(Histogram,slens,0:binstep:xmax))
    lens = DataFrame(height = [hheight.weights; sheight.weights],
                     experiment = repeat(["human","simulation"],
                                         inner=length(hheight.weights)),
                     length = repeat((0:binstep:xmax)[1:end-1].+binstep/2,2))

    hplot = plot(lens,x=:length,color=:experiment,y=:height,
                 ygroup=:experiment,
                 Guide.colorkey(pos=[0.65*Gadfly.w,-0.3*Gadfly.h]),
                 Scale.color_discrete_manual("darkgray","lightgray"),
                 Geom.subplot_grid(Guide.xlabel("log-z-scored length",
                                                orientation=:horizontal),
                                   Guide.ylabel("density",orientation=:vertical),
                                   Coord.cartesian(xmin=0,xmax=xmax),
                                   Geom.bar,free_x_axis=true),
                 Theme(discrete_highlight_color=x->"black",
                       bar_highlight=x->"black"))
  end
end

function plot_fit(df,params;showci=true,
                  alpha=0.05,resample=400,kwds...)
  hstack(plot_stream(df,params;kwds...),
         plot_lengths_qq(df,params;kwds...,showci=showci,
                         alpha=alpha,resample=resample))
end

function plot_responses((len,value),args...)
  dfl = DataFrame(value=value,time=cumsum(len));
  dfl.lagtime = lag(dfl.time,default=0.0)
  dfl[!,:ymin] .= -1
  dfl[!,:ymax] .= 1

  plot(dfl,xmax=:time,ymin=:ymin,xmin=:lagtime,ymax=:ymax,color=:value,
       Scale.color_discrete_manual("black","lightgray"),
       Geom.rect,args...)
end

function plot_fitmask(data,params,settings;Δf=6,simulation=1,start_time=0s,
                      stop_time=20s,kwds...)
  fit = plot_fit(data,params;kwds...)
  mask = plot_mask(data,params,settings;Δf=Δf,simulation=simulation,
                   start_time=start_time,stop_time=stop_time,kwds...)
  vstack(fit,mask)
end

function plot_mask(df,params,settings;Δf=6,simulation=1,start_time=0s,
                   stop_time=20s,kwds...)
  selection = select_params(params;Δf=Δf,kwds...)
  mask = select_mask(df,params,settings;Δf=Δf,simulation=simulation,
                     start_time=start_time,stop_time=stop_time,kwds...)
  input = audiospect_stimulus(params[selection,:],settings)
  input = input[start_time .. stop_time]
  band = plot_responses(percept_lengths(mask,input,settings))

  spect=plot(asplotable(mask,quantize=(200,128))[1],
             x=:time,y=:logfreq,color=:value,Geom.rectbin,
             Coord.cartesian(xmin=ustrip(start_time),xmax=ustrip(stop_time)),
             Scale.color_continuous(colormap=Scale.lab_gradient("white","black")))

  vstack(band,spect)
end
