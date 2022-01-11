using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
include("MPO.jl")
include("variational_first_excited_state_MPS_algorithm.jl")

# ----------------------------------------------------------------------------------------------------------------------------------

# Generate data for extrapolation to continuum - # N_list = [40, 60, 80, 100], D_list = [20, 40, 60, 80, 100], x_list = [30.0, 40.0, 50.0, 60.0, 80.0] # N/sqrt(x) > 5, mg = 0.125

# accuracy = 10^(-8)
# lambda = 100.0
# l_0 = 0.0
# max_sweep_number = 100

# N = parse(Int, ARGS[1])
# D = parse(Int, ARGS[2])
# D_previous = parse(Int, ARGS[3])
# x = parse(Float64, ARGS[4])
# mg = parse(Float64, ARGS[5])

# generate_Schwinger_data(mg, x, N, D, accuracy, lambda, l_0, max_sweep_number, D_previous)

# ----------------------------------------------------------------------------------------------------------------------------------

# Generate data for entanglement entropy vs mass plot 

generate_entropy_data()

# ----------------------------------------------------------------------------------------------------------------------------------
