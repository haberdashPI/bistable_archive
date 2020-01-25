using Test
using CorticalSpectralTemporalResponses
using AuditoryBistabilityLE
using SampledSignals
using TOML
using LinearAlgebra
using Statistics
using JLD2
using FileIO
using DataFrames
using Unitful

N = 50

pstream(p) = sum(p.counts[1][p.counts[2]]) / sum(p.counts[1])

@testset "Stremaing Bistability" begin
  @testset "Basic Streaming" begin
    params = load("object_test_params.jld2")["params"]

    # TODO: get this function running

    settings = TOML.parsefile("test_settings.toml")
    settings["stimulus"]["repeats"] = 10

    @test bistable_model(params, settings, interactive=true,
                   progressbar=false) != false
    # uncomment to profile
    # @time bistable_model(100, params, "test_settings.toml", interactive=true,
    #                      progressbar=true)

    # TODO: test with the always-fused stimulus from Elhilali et al 2009
    # (not at immediate concern, but should be present for the final
    # mod)
    params[:Δf] = 12
    settings["stimulus"]["repeats"] = 50
    df12 = bistable_model(params, "test_settings.toml",
                          interactive=true, progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(params, "test_settings.toml", interactive=true,
                         progressbar=true)

    @test pstream(df12.percepts) > 0.9
    @test pstream(df3.percepts) < 0.1
  end

  @testset "Track-level Bistability" begin

    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :f_c_a => 0, :f_c_m => 0, :f_c_σ => 0,
      :s_c_a => 0, :s_c_m => 0, :s_c_σ => 0,
      :t_W_m_σ      => 15.0, #5.0
      :t_W_m_σ_t    => 7.0,   :t_W_m_σ_ϕ   => 7.0,
      :t_W_m_σ_N    => 3.0,   :t_W_m_c     => 6.0,
      :t_τ_m        => 350ms, :t_c_m       => 100,
      :t_τ_a        => 3s,    :t_c_a       => 6,
      :t_τ_σ        => 500ms, :t_c_σ       => 0.2
     )

    params[:Δf] = 12
    df12 = bistable_model(100, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 6
    df6 = bistable_model(100, params, "test_settings.toml", interactive=true,
                         progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(100, params, "test_settings.toml", interactive=true,
                         progressbar=true)


    # bare minimum tests that should pass if there's stimulus selective
    # bistability the results are somewhat random, so not all tests will pass
    # on all runs, but the should normally *mostly* pass.
    # once the parameters have been fit well, the goal will be to select
    # parameter values that consistently pass
    @test pstream(df12.percepts) > 0.8
    @test pstream(df3.percepts) < 0.2
    @test pstream(df12.percepts) > pstream(df6.percepts)
    @test pstream(df6.percepts) > pstream(df3.percepts)
  end

  @testset "Scale-level Bistability" begin

    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :f_c_a => 0, :f_c_m => 0, :f_c_σ => 0,
      :s_W_m_σ      => 15.0,
      :s_W_m_c      => 6.0,
      :s_τ_m        => 350ms, :s_c_m       => 0,
      :s_τ_a        => 3s,    :s_c_a       => 25,
      :s_τ_σ        => 500ms, :s_c_σ       => 0.2,
      :t_c_a => 0, :t_c_m => 0, :t_c_σ => 0
     )

    params[:Δf] = 12
    df12 = bistable_model(100, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 6
    df6 = bistable_model(100, params, "test_settings.toml", interactive=true,
                         progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(100, params, "test_settings.toml", interactive=true,
                         progressbar=true)


    # bare minimum tests that should pass if there's stimulus
    # selective bistability
    @test pstream(df12.percepts) > 0.8
    @test pstream(df3.percepts) < 0.2
    @test pstream(df12.percepts) > pstream(df6.percepts)
    @test pstream(df6.percepts) > pstream(df3.percepts)
  end


  @testset "Freq-level Bistability" begin
    params = Dict(
      :Δt         => 240ms, :Δf        => 3,
      :f          => 500Hz, :condition => :track,
      :τ_x        => 500ms, :c_x       => 3.0,
      :f_W_m_σ      => 5.6,
      :f_W_m_c      => 6.0,
      :f_τ_m        => 350ms, :f_c_m       => 32,
      :f_τ_a        => 3s,    :f_c_a       => 5,
      :f_τ_σ        => 500ms, :f_c_σ       => 0.2,
      :s_c_a => 0, :s_c_m => 0, :s_c_σ => 0,
      :t_c_a => 0, :t_c_m => 0, :t_c_σ => 0
     )

    params[:Δf] = 12
    df12 = bistable_model(50, params, "test_settings.toml", interactive=true,
                          progressbar=true)

    params[:Δf] = 6
    df6 = bistable_model(50, params, "test_settings.toml", interactive=true,
                         progressbar=true)

    params[:Δf] = 3
    df3 = bistable_model(50, params, "test_settings.toml", interactive=true,
                          progressbar=true)


    # bare minimum tests that should pass if there's stimulus
    # selective bistability
    @test pstream(df12.percepts) > 0.8
    @test pstream(df3.percepts) < 0.2
    @test pstream(df12.percepts) > pstream(df6.percepts)
    @test pstream(df6.percepts) > pstream(df3.percepts)
  end
end

