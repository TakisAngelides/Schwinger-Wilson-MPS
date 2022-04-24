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

# Generate data for entanglement entropy vs mass plot 

N = parse(Int, ARGS[1])
x = parse(Float64, ARGS[2])
mg = parse(Float64, ARGS[3])
D = parse(Int64, ARGS[4])
l_0 = parse(Float64, ARGS[5])
lambda = parse(Float64, ARGS[6])
accuracy = parse(Float64, ARGS[7])
max_sweep_number = parse(Int64, ARGS[8])

generate_entropy_data(mg, x, N, D, accuracy, lambda, l_0, max_sweep_number)

# ----------------------------------------------------------------------------------------------------------------------------------
