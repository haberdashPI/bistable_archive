include("base_setup.jl")
using RCall
R"""
library(dplyr)
library(cowplot)
library(ggplot2)

options(rcalljl_options=list(width=800,height=400))
"""

# import Cairo, Fontconfig
using Revise # (we have to call this in IJulia... seems like a bug)
revise()

include(joinpath(srcdir,"count_lengths.jl"))
include(joinpath(srcdir,"measures.jl"))
include(joinpath(srcdir,"plotting.jl"))


