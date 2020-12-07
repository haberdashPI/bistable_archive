using Pkg
Pkg.activate(joinpath(@__DIR__,".."))

using Unitful
using AxisArrays
using TOML
using LinearAlgebra
using Statistics
using Dates
using ProgressMeter
using Colors
using Formatting
using CSV
using DependentBootstrap
using FFTW

const srcdir = joinpath(pwd(),"..", "src")
const plotdir = joinpath(pwd(),"..","plots","paper")
const grantdir = joinpath(pwd(),"..","plots","grant")

