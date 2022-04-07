using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
include("MPO.jl")
include("variational_first_excited_state_MPS_algorithm.jl")

accuracy = 10^(-8)
lambda = 100.0
l_0 = 0.0
max_sweep_number = 100
N = parse(Int, ARGS[1])
x = parse(Float64, ARGS[2])
mg = parse(Float64, ARGS[3])
D = parse(Int64, ARGS[4])

mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)

spin_N = 2*N

path_to_mps = "/lustre/fs23/group/nic/tangelides/Schwinger Wilson MPS Data/N_$(spin_N)_x_$(x)_D_$(D)_l0_$(l_0)/mps_$(spin_N)_$(D)_$(mg)_$(x)_$(l_0).h5"

f = h5open(path_to_mps, "r")

mps_group = f["$(lambda)_$(l_0)_$(mg)_$(x)_$(spin_N)_$(D)"]

mps_ansatz = Vector{Array{ComplexF64}}(undef, spin_N)

for i in 1:spin_N

    mps_ansatz[i] = read(mps_group["mps_$(i)"])

end

_, _, _ = variational_ground_state_MPS_from_dead_for_saving(2*N, 2, D, mpo, accuracy, max_sweep_number, mps_ansatz)