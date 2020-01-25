module AuditoryBistabilityLE
using CorticalSpectralTemporalResponses
using DSP
using DataFrames
using Requires
using Statistics
using LinearAlgebra
using SparseArrays
using SpecialFunctions
using TOML
using MetaArrays
using Unitful

import CorticalSpectralTemporalResponses: Δt, Δf, times, freqs, scales, rates, asseconds, asHz

export adaptmi, drift, scale_weighting, ncomponents, nunits, CoherenceModel,
    fusion_ratio, object_SNR, mask, scene_object_ratio,
    object_SNR2, ab_match, mean_spect, mean_spect2, AdaptMI

using ProgressMeter
next!(x::Progress) = ProgressMeter.next!(x)
next!(x::Nothing) = nothing

# core model implementation
########################################
# implements inhibition, adaptation and noise
include("adaptmi.jl")
# implements inhibition interactions for each hierarchical level
include("cortmi.jl")
# utility functions for handling the windowing of a signal
include("windowing.jl")
# cohere is step 1 of the object-level analysis (temporal coherence)
include("cohere.jl")
# the parts of cohere specific to the NMF implementation
include("nmf.jl")
# general purpose functions for tracking (step 2 of the object-level analysis)
include("tracking.jl")
# function specific to using prior source models for object tracking
# (other approachs were examined)
include("prior_tracking.jl")
# function for considering multiple priors
# (i.e. the different β, κ parametesr from the paper)
include("multi_prior_tracking.jl")

# all the functions that use the model implementation to do stuff
########################################
# miscellaneous utility functions
include("util.jl")
# ABA stimulus generation
include("stim.jl")
# hueristic estimate of the number of sources present in an ABA_ stimulus
include("bandwidth_ratio.jl")
#  code to compute percept lengths using bandwidth ratio
include("lengths.jl")
# applies the bistable model to a specific set of conditions
include("biapply.jl")
# where everything comes together: all stages of the model are computed here
include("bimodel.jl")
# compress the model output so it takes up less space when stored
include("compress.jl")

# data visualization
include("plot_axes.jl")
function __init__()
  @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("rplots.jl")
end

end
