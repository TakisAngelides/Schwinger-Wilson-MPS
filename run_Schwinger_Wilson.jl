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

# D values to run: 20 40 60 80 100 110 120 140 160 180 200 250 300
# mg values to run: -0.10, -0.20, -0.30, -0.45, -0.47, -0.49, -0.51, -0.53, -0.55, -0.60, -0.65, -0.75
# extra mg values for N = 48, x = 10.0: -0.56, -0.57, -0.58, -0.59, -0.61, -0.62, -0.63, -0.64, -0.66, -0.67, -0.68, -0.69, -0.71, -0.72  
accuracy = 10^(-8)
lambda = 100.0
l_0 = 0.5
max_sweep_number = 100
N = parse(Int, ARGS[1])
x = parse(Float64, ARGS[2])
mg = parse(Float64, ARGS[3])
D = parse(Int64, ARGS[4])
choice = parse(Int64, ARGS[5])

if choice == 1
    generate_entropy_data(mg, x, N, D, accuracy, lambda, l_0, max_sweep_number)
elseif choice == 2
    mps_to_entropy_save_file(mg, x, 2*N, D)
elseif choice == 3
    mps_to_average_electric_field(mg, x, 2*N, D, l_0)
end

# ----------------------------------------------------------------------------------------------------------------------------------
