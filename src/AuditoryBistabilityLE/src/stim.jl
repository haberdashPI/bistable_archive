using Unitful
using SampledSignals

export aba_, ab, stimulus, audiospect_stimulus

const sr = 8000

silence(len) = SampleBuf(zeros(SampledSignals.inframes(Int,len,sr)),sr)
function tone(freq,len) 
  freq_ = SampledSignals.inHz(freq)
  len_ = SampledSignals.inseconds(len)
  SampleBuf(sin.(2π.*freq_.*range(0,stop=len_,step=1/sr)), sr)
end

function ramp(xs,len)
  n = SampledSignals.inframes(Int,len,samplerate(xs))
  t = range(0,stop=1,length=n)
  xs[1:n] .*= cos.(π.*(t.-1)./2)
  xs[end-n+1:end] .*= cos.(π.*t./2)
  xs
end

normpower(x) = x ./ sqrt.(mean(x.*x,dims=1))
amplify(ratio) = x -> x.*ratio

const stim_cache = Dict{Vector{Union{Number,String}},AbstractMatrix}()
function audiospect_stimulus(params,settings;cache=false)
  settings = read_settings(settings)
  params = read_params(params)
  if cache
    get!(stim_cache,[[params[:Δt],params[:Δf],params[:f]];
                             collect(values(settings.stimulus))]) do
      audiospect_stimulus_(params,settings)
    end
  else
    audiospect_stimulus_(params,settings)
  end
end

function audiospect_stimulus_(params,settings)
  stim = stimulus(params[:Δt],params[:f],params[:Δf];
                  settings.stimulus...) |> normpower |> amplify(-10dB)
  @info "Stimulus is $(maximum(domain(stim))) seconds long."
  as = Audiospect(;settings.freqs.analyze...)
  filt(as,stim,false)
end

function stimulus(total_len,freq,delta;repeats=10,tone_len_fraction=0.5,
                  pattern="ab",ramp_len=0ms)
  if pattern == "ab"
    ab(tone_len_fraction*total_len,(1-tone_len_fraction)*total_len,1,
       repeats,freq,delta;ramp_len=asseconds(ramp_len))
  elseif pattern == "aba_"
    aba_(tone_len_fraction*total_len,(1-tone_len_fraction)*total_len,
         repeats,freq,delta;ramp_len=asseconds(ramp_len))
  else
    error("Unexpected stimulus pattern '$pattern'")
  end
end

function aba_(tone_len,spacing_len,repeats,freq,delta;ramp_len=0ms)
  a_freq=freq
  b_freq=freq*(2^(delta/12))
  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  a = ramp_len > 0ms ? ramp(a,25ms) : a
  b = ramp_len > 0ms ? ramp(b,25ms) : b

  aba_ = [a;silence(spacing_len);
          b;silence(spacing_len);
          a;silence(2spacing_len+tone_len)]

  aba_seq = silence(0s)
  for i = 1:repeats
    aba_seq = [aba_seq; aba_]
  end

  aba_seq
end

function ab(tone_len,spacing_len,offset_ratio,repeats,freq,delta;ramp_len=0ms)
  @assert 0 <= offset_ratio <= 1

  a_freq=freq
  b_freq=freq*(2^(delta/12))
  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  a = ramp_len > 0ms ? ramp(a,25ms) : a
  b = ramp_len > 0ms ? ramp(b,25ms) : b

  offset = silence(offset_ratio*(tone_len+spacing_len))
  ab = [a;offset].+[offset;b]
  ab_len = 2(tone_len+spacing_len)
  ab = [ab; silence(ab_len - duration(ab)*s)]

  ab_seq = silence(0s)
  for i = 1:repeats
    ab_seq = [ab_seq; ab]
  end

  ab_seq
end
