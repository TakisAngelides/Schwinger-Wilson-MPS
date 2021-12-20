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

# Checking that the spectrum of the exact Schwinger Hamiltonian agrees with the matrix formed by its MPO

# N = 2
# x = 1.0
# m_g_ratio = 0.5
# l_0 = 0.0 # this is theta/2pi
# lambda = 0.0

# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)

# matrix = mpo_to_matrix(mpo)

# matrix_h = get_Schwinger_hamiltonian_matrix(N, l_0, x, lambda, m_g_ratio)

# display(norm(eigvals(matrix_h)-eigvals(matrix)))

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking that the minimum energy from the variational ground state search agrees with the minimum energy from exact diagonalization

# N = 2
# x = 1.0
# m_g_ratio = 0.4
# l_0 = 0.5 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 20
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)
# println("Minimum energy from variational ground state search: ", E_0)
# matrix = mpo_to_matrix(mpo)
# display(eigvals(matrix))
# println("Minimum energy from exact diagonalization: ", minimum(eigvals(matrix)))

# ----------------------------------------------------------------------------------------------------------------------------------

# Check that the penalty term enforces total charge to 0 and checking the local charge MPO

# N = 4
# x = 1.0
# m_g_ratio = 0.5
# l_0 = 0.25 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 10
# d = 2
# D = 10

# penalty_mpo = get_penalty_term_MPO(N, lambda)
# penalty_mpo_matrix = mpo_to_matrix(penalty_mpo)
# penalty_matrix_exact = get_penalty_term_matrix(N, lambda)
# display(eigvals(penalty_mpo_matrix))
# display(eigvals(penalty_matrix_exact))

# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)

# mpo_penalty = get_penalty_term_MPO(N, lambda)
# mps_after_penalty_mpo = act_mpo_on_mps(mpo_penalty, mps_ground)

# println(get_mpo_expectation_value(2*N, mps_ground, mpo_penalty))
# println(get_spin_half_expectation_value(2*N, mps_ground, "z")) # the total charge operator is sum -g/2 sigma_z

# charge_list = []

# for n in 2:2:2*N

#     charge_mpo = get_local_charge_MPO(N, n)
#     mps_right = act_mpo_on_mps(charge_mpo, mps_ground)
#     append!(charge_list, inner_product_MPS(mps_ground, mps_right))

# end

# display(charge_list)
# println(sum(charge_list))

# l_field = get_electric_field_configuration(N, l_0, mps_ground)

# # This checks Gauss' law: l_field[i] - l_field[i-1] = Q[i]

# display(isapprox(l_field[2]-l_field[1], charge_list[2]))
# display(l_field[2]-l_field[1])
# display(charge_list[2])

# display(sum(get_electric_field_configuration(N, l_0, mps_ground))/N)

# ----------------------------------------------------------------------------------------------------------------------------------

# Check that the first excited energy is the same as from exact diagonalization

# N = 4
# x = 1.0
# m_g_ratio = 0.125
# l_0 = 0.0
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 20
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)
# E_1, mps_first, sn_first = variational_first_excited_MPS(2*N, d, D, mpo, acc, max_sweeps, mps_ground, E_0)
# println("First excited energy from variational search: ", E_1)
# matrix = mpo_to_matrix(mpo)
# println("First excited energy from exact diagonalization: ", eigvals(matrix))

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking memory allocations

# # The command which you write in a terminal not in REPL to generate the variational_MPS_algorithm.jl.mem file is:

# # julia --track-allocation=user variational_MPS_algorithm.jl

# # Then you run the variational_MPS_algorithm.jl again in the terminal with the command:

# # julia variational_MPS_algorithm.jl 

# # and then open the .mem file which will contain the number of memory allocations (in bytes?) per line of code.

# function wrapper() # so as to not misallocate and focus on the function we want to probe
#     mpo = get_Schwinger_Wilson_MPO(80, 0.0, 1.0, 100.0, 0.125) # force compilation
#     Profile.clear_malloc_data() # clear allocation
#     mpo = get_Schwinger_Wilson_MPO(80, 0.0, 1.0, 100.0, 0.125) # run again without compilation
# end

# wrapper()

# ----------------------------------------------------------------------------------------------------------------------------------

# Check the contents of an h5 file

# h5open("mps_60_20_0.125_30.0.h5", "r") do fid

#     g = fid["100.0_0.0_0.125_30.0_60_20"]

#     println(length(read(g)))
#     # for t in g
#     #     println(size(read(t)))
#     # end
#     # display(size(read(g["mps_120"])))
#     a = read(g["mps_1"])
#     display(a)
#     println(Base.summarysize(a)) # to print the memory needed for the variable a in bytes
# end

# ----------------------------------------------------------------------------------------------------------------------------------

# # https://github.com/kuehnste/mps.jl/blob/main/test/runtests.jl

# @testset "MPS Testing" begin

#     @test 

# end

# ----------------------------------------------------------------------------------------------------------------------------------

