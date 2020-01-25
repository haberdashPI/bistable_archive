export map_components
using Base.Iterators: product

@with_kw struct MultiPriorTracking{C} <: Tracking
  cohere::C
  time_constants::Array{typeof(1.0s)}
  time_constant_bias::Vector{Float64}
  source_priors::AxisArray
  freq_prior
  max_sources::Int = 4
  min_norm::Float64 = Inf
  normalize::Bool = false
end
function Tracking(C,::Val{:multi_prior};
                  time_constants=[4s],
                  time_constant_bias=zeros(length(time_constants)),
                  source_prior_sds,
                  source_prior_strengths=[1.0s],
                  source_prior_sd_bias=zeros(length(source_prior_sds)),
                  source_prior_strength_bias=
                    zeros(length(source_prior_strengths)),
                  freq_ridge=0.0, scale_ridge=0.0,
                  ridge_threshold=0.05,
                  freq_prior_N = 2, freq_prior_bias = 0,
                  normalize=false,min_norm=Inf,
                  params...)
  checknorm(normalize,min_norm)
  source_prior_Ns = asseconds.(source_prior_strengths)./Δt(C)
  pr = product(zip(source_prior_sds, source_prior_sd_bias),
               zip(source_prior_Ns, source_prior_strength_bias))

  source_priors = if iszero(freq_ridge) && iszero(scale_ridge)
    vals = [(isonorm(sd, N, (size(C,2),size(C,3))), sd_b + N_b)
                            for ((sd,sd_b),(N,N_b)) in pr]
    prior_params = map(x -> map(a -> a[1],x),pr)
    AxisArray(vec(vals), Axis{:params}(vec(prior_params)))
  else
    vals = [(ridgenorm(sd, N, (size(C,2),size(C,3)),
                       freq=freq_ridge,scale=scale_ridge,
                       threshold=ridge_threshold), sd_b + N_b)
            for ((sd,sd_b),(N,N_b)) in pr]
    prior_params = map(x -> map(a -> a[1],x),pr)
    AxisArray(vec(vals), Axis{:params}(vec(prior_params)))
  end

  freq_prior = freqprior(freq_prior_bias,freq_prior_N)
  MultiPriorTracking(;cohere=getmeta(C), source_priors=source_priors,
                     freq_prior=freq_prior,
                     time_constants=asseconds.(time_constants),
                     time_constant_bias=time_constant_bias,
                     normalize=normalize,min_norm=min_norm,
                     params...)
end

function expand_params(params::MultiPriorTracking)
  pr = params.source_priors

  AxisArray([(PriorTracking(params.cohere,tc,prior,params.freq_prior,
                            params.max_sources,params.min_norm), tb + pb)
             for (tc,tb) in zip(params.time_constants,
                                params.time_constant_bias)
             for (prior,pb) in pr],
            Axis{:params}([(tc,sd,N) for tc in params.time_constants
                           for (sd,N) in axisvalues(params.source_priors)[1]]))
end

function nitr(C::Coherence,params::MultiPriorTracking)
  ntimes(C) * length(params.time_constants) * length(params.source_priors)
end

function track(C::Coherence,params::MultiPriorTracking,progressbar=true,
               progress=track_progress(progressbar,nitr(C,params),"multi-prior"))

  C_ = prepare_coherence(C,params.min_norm)
  all_params = expand_params(params)
  S = Array{SourceTracking}(undef,size(all_params,1))
  lp = Array{Array{Float64}}(undef,size(all_params,1))

  # in my tests, using @Threads does not seem to help hereo This might be
  # related to a now fixed bug in Julia
  # (https://github.com/JuliaLang/julia/issues/17395) or it could because the
  # operation is memory bound; in either case, I haven't gone back to test this
  # and see if threading works any better now
  #=@Threads.threads=# for (i,(p,bias)) in collect(enumerate(all_params))
    S[i], lp[i] = track(C_,p,true,progress)
    lp[i] .+= bias
  end

  (AxisArray(S, AxisArrays.axes(all_params,1)),
   AxisArray(hcat(lp...), AxisArrays.axes(C,1), AxisArrays.axes(all_params,1)))
end

function map_components(fn,tracks::AxisArray{<:SourceTracking},
                        tracks_lp::AxisArray{<:Float64};
                        window=500ms,step=250ms)
  windows = windowing(tracks[1],length=window,step=step)

  result = map(enumerate(windows)) do (i,ixs)
    best_track = argmax(dropdims(mean(Array(tracks_lp[ixs,:]),dims=1),dims=1))
    fn(tracks[best_track][Axis{:time}(ixs)])
  end

  AxisArray(result,AxisArrays.axes(windows,Axis{:time}))
end

function mask(sp::CorticalSpectralTemporalResponses.AuditorySpectrogram,
              tracks::AxisArray{<:SourceTracking},
              tracks_lp::AxisArray{<:Float64},
              settings;progressbar=false,kwds...)
  scales = settings["scales"]["values"] .* cycoct
  freql,freqh = settings["rates"]["freq_limits"] .* Hz

  cr = cortical(sp,scales=scales,progressbar=progressbar)
  cr = cr[:,:,freql .. freqh]

  mask(cr,tracks,tracks_lp;progressbar=progressbar,kwds...)
end

function mask(cr::MetaUnion{AxisArray},
              tracks::AxisArray{<:SourceTracking},
              tracks_lp::AxisArray{<:Float64},
              order=1;window=500ms,mask_wait=0,
              delta=250ms,progressbar=false)

  @assert axisdim(cr,Axis{:time}) == 1
  @assert axisdim(cr,Axis{:scale}) == 2
  @assert axisdim(cr,Axis{:freq}) == 3
  @assert size(cr)[2:3] == size(tracks[1])[1:2] "Dimension mismatch"

  windows = windowing(tracks[1],length=window,step=delta)

  progress = progressbar ? Progress(length(windows),"Masking: ") : nothing
  mask_helper(cr,tracks,tracks_lp,order,windows,mask_wait,progress)
end

function mask_helper(cr,tracks,tracks_lp,order,windows,mask_wait,progress)
  y = zero(cr)
  norm = similar(cr,real(eltype(cr)))
  norm .= zero(eltype(norm))

  cohere_windows =
    collect(windowing(cr,tracks[1].cohere))

  for (i,ixs) = enumerate(windows)
    best_track = argmax(dropdims(mean(Array(tracks_lp[ixs,:]),dims=1),dims=1))
    components = view(tracks[best_track],:,:,:,ixs)
    sorting = sortperm(component_means(components),rev=true)
    c = sorting[order]

    for (ti,t) in enumerate(ixs)
      if i <= mask_wait
        y[cohere_windows[t],:,:] .+= 1
        norm[cohere_windows[t],:,:] .+= 1
      else
        resh = reshape(view(components,:,:,c,ti),1,size(y,2),size(y,3))
        y[cohere_windows[t],:,:] .+= resh
        norm[cohere_windows[t],:,:] .+= sqrt.(mean(x -> x^2,resh))
      end
    end

    next!(progress)
  end
  y ./= max.(1,norm)
  y ./= maximum(abs,y)
  y .*= cr

  y
end
