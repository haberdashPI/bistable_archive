using Pkg
Pkg.activate(joinpath(@__DIR__,".."))

using Unitful
using AxisArrays
using TOML
using DataFramesMeta
using LinearAlgebra
using Statistics
using Dates
using ProgressMeter
using Colors
using Formatting
using CSV
using DependentBootstrap
using FFTW
using Query

const srcdir = @__DIR__
const plotdir = joinpath(@__DIR__,"..","plots","paper")
const grantdir = joinpath(@__DIR__,"..","plots","grant")

