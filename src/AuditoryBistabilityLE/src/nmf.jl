using NMF
using Parameters

@with_kw struct CoherenceNMF
  normalize::Bool = true
  maxiter::Int = 2000
  tol::Float64 = 1e-4
  warn_converge::Bool = true
  reuse_solution::Bool = false
end
CoherenceMethod(::Type{Val{:nmf}},params) = CoherenceNMF(;params...)

Base.eltype(::CoherenceNMF,::Type) = Float64
Base.eltype(::CoherenceNMF,::AbstractArray) = Float64

function with_method(foreach_window,method::CoherenceNMF,K)
  convergence_count,total_count = 0,0
  Winit,Hinit = fill(0.0,(0,0)),fill(0.0,(0,0))
  method_ = NMF.ALSPGrad{Float64}(tol=method.tol,
                                 maxiter=method.maxiter)

  foreach_window() do window
    x_t = reshape(window,size(window,1)*size(window,2),:)
    k = min(size(x_t,1),K)
    nmf_input = abs.(x_t) .+ method.tol/2

    if isempty(Winit) || !method.reuse_solution
      Winit,Hinit = NMF.nndsvd(nmf_input,k)
    end

    solution = NMF.solve!(method_,nmf_input,Winit,Hinit)
    total_count += 1

    # solution = nnmf(abs.(x_t) .+ method.tol/2,k,init=:nndsvd,
    #                 tol=method.tol,
    #                 maxiter=method.maxiter)

    convergence_count += !solution.converged

    y = reshape(solution.H,:,size(window)[3:end]...)
    permutedims(y,[2:ndims(window)-1;1])
  end

  if convergence_count > 0 && method.warn_converge
    percent = round(100convergence_count / total_count,digits=2)
    if percent > 0
      @info("$percent% of frames (n = $convergence_count) failed to fully "*
            "converge to a solution.")
    else
      @info("<0.01% of frames (n = $convergence_count) failed to fully "*
            "converge to a solution.")
    end
  end
end
