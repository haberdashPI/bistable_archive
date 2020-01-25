using RCall
include("../../models/little/count_lengths/count_lengths.jl")

signal = AxisArray(ones(2000,2),Axis{:time}(linspace(0,2,2000)*s))
drifted = drift(signal,τ_σ=0.5s,c_σ=0.5)

R"""
pdf("noise.pdf")
plot($(Array(drifted[:,1])),type='l')
dev.off()
"""

