using StatsBase: fit, Histogram, weights
using Random

rms(x) = sqrt(mean(x.^2))
meansqr(x) = mean(x.^2)

function select_params(params;kwds...)
  condition = trues(size(params,1))
  for (var,val) in pairs(kwds)
    condition .&= ustrip.(abs.(params[:,var] .- val)) .<= 1e-2
  end
  params[condition,:pindex]
end

function select_data(df,params;kwds...)
  sel = select_params(params;kwds...)
  dfsel = df[in.(df.pindex,Ref(sel)),:]
  found = params[unique(dfsel.pindex),:]
  if size(found,1) < length(sel)
    kwdbuf = IOBuffer()
    show(kwdbuf,kwds.data)
    @warn("Expected $(length(sel)) parameter entries. "*
          "\nInstead, only found entires: "*string(found)*
          "\nKeyword selection: "*String(take!(kwdbuf)))
  end
  dfsel.st = params.Δf[indexin(dfsel.pindex,params.pindex)]
  dfsel, params[sel,:]
end

function model_data(df,params;bound=true,kwds...)
  df,params = select_data(df,params;kwds...)
  (stream=stream_summary(df,params,bound=bound), lengths=length_summary(df,params))
end

error_ratio(a,b=human_error()) = (a.stream / b.stream + a.lengths / b.lengths)/2
function model_error(df::DataFrame,params::DataFrame;bound=true,kwds...)
  mdata = model_data(df,params;bound=bound,kwds...)
  model_error(mdata,data_summarize(human_data()))
end

function model_error(data::NamedTuple,mean::NamedTuple)
  str_error = mean_bysid(data.stream) do str
    stream_rms(str,mean.stream)^2
  end |> sqrt
  len_error = ksstat(data.lengths,mean.lengths.lengths)

  (stream=str_error,lengths=len_error)
end

function model_stream_stats(data::DataFrame,params::DataFrame;kwds...)
  mdata = model_data(data,params;kwds...)
  model_stream_stats(mdata,data_summarize(human_data()))
end

function model_stream_stats(data::NamedTuple,mean::NamedTuple)
  shift = mean_bysid(data.stream) do str
    stream_shift(str,mean.stream)
  end
  spread = mean_bysid(data.stream) do str
    stream_spread(str,mean.stream)
  end
  (shift=shift,spread=spread)
end

function mean_bysid(fn,data)
  total = 0.0
  n = 0
  for g in groupby(data,:sid)
    total += fn(g)
    n += 1
  end
  total / n
end

function stream_shift(data,mean)
  diffs = map(1:size(mean,1)) do row
    i = findfirst(isequal(mean.st[row]),data.st)
    if i isa Nothing
      missing
    else
      data[i,:streaming] - mean[row,:streaming]
    end
  end
  sum(coalesce.(diffs,0.0))
end

const st_spread = Dict(3 => -1,6 => 0,12 => 1)
function stream_spread(data,mean)
  diffs = map(1:size(mean,1)) do row
    i = findfirst(isequal(mean.st[row]),data.st)
    if i isa Nothing
      missing
    else
      st_spread[mean.st[row]]*(data[i,:streaming] - mean[row,:streaming])
    end
  end
  sum(coalesce.(diffs,0.0))
end

function stream_rms(data,mean)
  diffs = map(1:size(mean,1)) do row
    i = findfirst(isequal(mean.st[row]),data.st)
    if i isa Nothing
      missing
    else
      data[i,:streaming] - mean[row,:streaming]
    end
  end
  if all(ismissing.(diffs))
    return missing
  else
    rms(coalesce.(diffs,0.0))
  end
end

function ksstat(x,y)
  nx = length(x)
  ny = length(y)
  order = sortperm([x; y])
  diffs = cumsum(map(i -> i <= nx ? 1/nx : -1/ny,order))
  maximum(abs,diffs)
end

function findci(x;kwds...)
  if length(x) > 3
    if any(xi != x[1] for xi in x)
      cint = dbootconf(collect(skipmissing(x));kwds...)
      (mean=mean(skipmissing(x)), lowerc=cint[1], upperc=cint[2])
    else
      (mean=x[1], lowerc=x[1], upperc=x[1])
    end
  else
    (mean = mean(x), lowerc = minimum(x), upperc = maximum(x))
  end
end

function stream_summary(data,params;bound=true,bound_threshold=0.8)
  result = by(data,[:st,:created]) do g
    DataFrame(streaming = streamprop(g.percepts,g.length,bound=bound,
                                     threshold=bound_threshold))
  end
  if all(ismissing,result.streaming)
    result = by(data,[:st,:created]) do g
      DataFrame(streaming = streamprop(g.percepts,g.length,bound=false))
    end
  end
  insertcols!(result, 1, sid = fill(0, size(result, 1)))
  for g in groupby(result,:st)
    g.sid .= 1:size(g,1)
  end
  sort!(result,(:sid,:st))

  result
end

function length_summary(data,params;Δf=6,bound_threshold=0.8)
  pindex = params.pindex[params.Δf .== Δf]
  mapreduce(vcat,groupby(data[data.pindex .== pindex,:],[:st,:created])) do g
    handlebound(ixs -> g.length[ixs],g.length,threshold=bound_threshold)
  end
  # data[data.pindex .== pindex,:length]
end

function mean_streaming(df;findci=false)
  by(df,:st) do st
    if findci
      DataFrame([findci(st.streaming)])
    else
      DataFrame(streaming=mean(st.streaming))
    end
  end
end

function human_data(;resample=nothing)
  (stream=human_stream_data(),lengths=human_length_data())
end

function data_summarize(data)
  stream = by(data.stream,:st) do g
    DataFrame(streaming = mean(skipmissing(g.streaming)))
  end
  lengths = data.lengths

  (stream=stream,lengths=lengths)
end

function human_error(;kwds...)
  str, len = human_error_by_sid(;kwds...)
  (stream=mean(skipmissing(str.x1)), lengths=mean(len.x1))
end

function human_error_by_sid()
  means,meanl = data_summarize(human_data())
  stream = human_stream_data()
  lengths = human_length_data()

  str_error = map(groupby(stream,[:sid,:experiment])) do str
    stream_rms(str,means)
  end |> combine
  len_error = combine(groupby(lengths,[:sid])) do len
    ksstat(len.lengths,meanl.lengths)
  end

  str_error, len_error
end

asnum(x::Real) = Float64(x)
asnum(x::String) = x == "NA" ? missing : parse(Float64,x)
function human_stream_data()
  df1 = CSV.read(joinpath(@__DIR__,"..","analysis","yerkes","stream_prop.csv"))
  df2 = CSV.read(joinpath(@__DIR__,"..","analysis","context","stream_prop.csv"))
  insertcols!(df2, 1, experiment = fill("3", size(df2, 1)))
  df = vcat(df1,df2)
  df.sid = string.(df.sid)
  df.streaming = asnum.(df.streaming)
  sort!(df,(:sid,:st))
  df
end

const N_for_pressnitzer_hupe_2006 = 23
const pressnitzer_hupe_binsize = 1/6

function human_length_data()
  df = CSV.read(joinpath(@__DIR__,"..","data","Higgins_et_al_unpublished","phasedurations.csv"))
  rename!(df,:phase => :lengths)
  df.lengths = asnum.(df.lengths)/10000 # from micorseconds to seconds
  rename!(df,:subject => :sid)
  df[df.code .== 140,:]
end

function filter_minlength(len,sid,x)
  x̃ = filter(x -> x ≥ len,x)
  if (1-length(x̃)/length(x)) > 0.05
    percent = (round(100(1-length(x̃)/length(x)),digits=1))
    @warn "Eliminating %$percent of data for $sid."
  end
  x̃
end

# THOUGHTS: in P&H 2006, the mean normalization is on a per-individual basis.
# This might incline us to use a mean on a per-simulation basis, but I think it
# makes more sense to do across all runs, because the simulation represents
# repeated measurements from the same "individual"
function norm_bymean(x;minlength = 0.1,sid = "UNKNOWN")
  x ./= mean(x)
end

function norm_bylogz(x;minlength = 0.1,sid = "UNKNOWN")
  x .-= mean(x)
  s = std(x)
  !iszero(s) ? exp.(x./s) : exp.(x)
end

function handlebound(fn,seconds;bound=true,threshold=0.8)
    if length(seconds) < 2
      fn(1:length(seconds))
    end

    if bound && (sum(seconds[2:end]) > threshold*sum(seconds))
      fn(2:length(seconds))
    else
      fn(1:length(seconds))
    end
end

function update_buildup!(value,time,run)
  j = 1
  ts = cumsum(run.length)
  for (i,t) in enumerate(time)
    while j <= Base.length(ts) && t > ts[j]
      j += 1
    end
    j <= Base.length(ts) || break
    value[i] += run.response[j]-1
  end
  value
end

function buildup_image(buildup_data;delta,length)
  len = length; length = Base.length

  runs = groupby(buildup_data,:run)
  times = range(0,len,step=delta)
  buildup = DataFrame(time=repeat(times,outer=length(runs)),
                      run=repeat(1:length(runs),inner=length(times)))
  buildup[:,:value] .= 0.0

  for (i,run) in enumerate(runs)
    build_run = view(buildup,buildup.run .== i,:)
    update_buildup!(build_run.value,build_run.time,run)
  end

  buildup
end

function buildup_mean(buildup_data;delta,length)
  buildup = DataFrame(time=range(0,length,step=delta))
  buildup[:,:value] .= 0.0
  runs = groupby(buildup_data,:run)
  for run in runs
    update_buildup!(buildup.value,buildup.time,run)
  end
  buildup.value ./= Base.length(runs)
  buildup
end

# sanity check
#=
switches = clamp.(randn(1000).*0.2 .+ 0.5,0,1)
plot(x=switches,Geom.density)

df = DataFrame(
  response = repeat([1,2],outer=1000),
  length = mapreduce(x -> [x,1-x],vcat,switches),
  run = repeat(1:1000,inner=2)
)
meandf = buildup_mean(df,delta=0.1,length=1.0)
plot(meandf,x=:time,y=:value,Geom.line)
=#

function streamprop(percepts,seconds;kwds...)
  handlebound(seconds;kwds...) do range
    sum(seconds[range][percepts[range] .== 2]) / sum(seconds[range])
  end
end

function stim_per_second(seconds;kwds...)
  handlebound(seconds;kwds...) do range
    length(range) / sum(seconds[range])
  end
end

