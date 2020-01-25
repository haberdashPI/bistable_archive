using Parameters
using Unitful: s, ms, ustrip

# question: should we somehow limit the unit response???
# e.g. with:
sig(x) = 1/(1+exp(-10(x-0.5)))

@with_kw struct AdaptMI{S,I}
  c_x::Float64 = 1.0
  τ_x::typeof(1.0s) = 300ms

  shape_y::S = identity

  c_a::Float64 = 5
  τ_a::typeof(1.0s) = 1.5s

  c_m::Float64 = 10
  τ_m::typeof(1.0s) = 50ms
  W_m::I = identity

  τ_σ::typeof(1.0s) = 100ms
  c_σ::Float64 = 0.3
end

################################################################################
# generic adaptation and mutual-inhibition operation
adaptmi(x;progressbar=true,kw...) = adaptmi(x,AdaptMI(;kw...),progressbar)
function adaptmi(x,params::AdaptMI,progressbar=true)
  @assert :time ∈ axisnames(x)
  time = Axis{:time}

  y = similar(x)

  shape_y = params.shape_y
  τ_x = params.τ_x; c_x = params.c_x
  τ_a = params.τ_a; c_a = params.c_a
  τ_m = params.τ_m; c_m = params.c_m; W_m = params.W_m

  # the (a)daptation (m)utual inhibition 
  a = similar(y) # a = adaptation
  m = similar(y) # m = mutual inhibition
  x̃ = similar(y) # (this is the smoothed version of x)

  # temporary, single time slice variables
  x̃_t = zero(y[time(1)])
  y_t = copy(x̃_t)
  a_t = copy(x̃_t)
  m_t = copy(x̃_t)

  dt_x = eltype(a_t)(Δt(x) / τ_x)
  dt_a = eltype(a_t)(Δt(x) / τ_a)
  dt_m = eltype(a_t)(Δt(x) / τ_m)

  progress = progressbar ? Progress(desc="Adapt/Inhibit: ",ntimes(y)) : nothing
  for ti in Base.axes(times(y),1)
    t = time(ti)

    @. begin # across all indices...
      # smooth input
      x̃_t += (c_x*x[t] - x̃_t)*dt_x

      # apply adaptation and inhibition
      y_t = shape_y((1 - c_a*a_t)*x̃_t - c_m*m_t)

      # update adaptation and inhibition
      a_t += (y_t - a_t)*dt_a
      m_t += ($W_m(y_t) - m_t)*dt_m
    end
    x̃[t] = x̃_t
    y[t] = y_t
    a[t] = a_t
    m[t] = m_t

    next!(progress)
  end
  y,a,m,x̃
end

################################################################################
# generic drifting noise function

drift(x,along_axes...;progressbar=true,kw...) =
  drift(x,AdaptMI(;kw...),along_axes,progressbar)
function drift(x,params::AdaptMI,along_axes=typeof.(AxisArayys.axes(x)),
    progressbar=true)

  τ_σ, c_σ = params.τ_σ, log(1+params.c_σ)
  time = Axis{:time}

  y = similar(x)
  σ_t = zero(x[time(1),(ax(1) for ax in along_axes)...])
  dims = size(σ_t)
  progress = progressbar ? Progress(desc="Drift: ",ntimes(y)) : nothing
  for t in Base.axes(times(y),1)
    ya_t = x[time(t)]

    σ_t .+= -σ_t.*(Δt(x)/τ_σ) .+ randn(dims).*(c_σ*sqrt(2Δt(x)/τ_σ))
    ya_t .*= exp.(σ_t)

    y[time(t)] = ya_t

    next!(progress)
  end
  y
end
