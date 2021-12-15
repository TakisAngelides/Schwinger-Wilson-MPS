using Profile
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

# Generate data for extrapolation to continuum - # N_list = [40, 60, 80, 100], D_list = [20, 40, 60, 80, 100]

# Generated data for N = 20, D = 40

# Generating data for N = 40, D = 20, 40, 60, 80, 100

accuracy = 10^(-10)
lambda = 100.0
l_0 = 0.0
max_sweep_number = 200
m_over_g_list = [0.125]
x_list = [30.0, 40.0, 50.0, 60.0, 80.0] # N/sqrt(x) > 5

N = parse(Int, ARGS[1])
D = parse(Int, ARGS[2])

generate_Schwinger_data(m_over_g_list, x_list, N, D, accuracy, lambda, l_0, max_sweep_number)

# ----------------------------------------------------------------------------------------------------------------------------------