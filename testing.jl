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

# Checking the center_orthogonalize! function

# mps = initialize_MPS(8, 2, 4)
# idx = 4
# println(inner_product_MPS(mps, mps))
# S = center_orthogonalize!(mps, idx)
# println(norm_center_orthogonal(mps, S, idx))

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking the entanglement entropy functions

# mps = get_GHZ_mps(8)
# display(quantum_state_coefficients(mps, 8))

# mpo = get_Schwinger_Wilson_MPO(8, 0.0, 1.0, 100.0, 0.33)
# acc = 10^(-10)
# ms = 30
# E_0, mps, sn = variational_ground_state_MPS(16, 2, 4, mpo, acc, ms)


# println("norm = ", inner_product_MPS(mps, mps))
# println("Trial = ", entanglement_entropy(mps, idx))
# println("Inefficient = ", entanglement_entropy_inefficient(mps, idx+1))
# println("EE = ", entanglement_entropy_old(mps, idx))
# S = center_orthogonalize!(mps, 4)
# evals = S.^2
# S_ortho = 0
# for eval in evals
#     if eval > 10^(-12)
#         global S_ortho += -eval*log2(eval)
#     end
# end
# println("S ortho = ", S_ortho)


# ----------------------------------------------------------------------------------------------------------------------------------

# Testing the GHZ state

# N = 6
# mps = get_GHZ_mps(N)
# mps = initialize_MPS(N, 2, 5)
# println(inner_product_MPS(mps, mps))
# println("EEI = ", entanglement_entropy_inefficient(mps, Int(N/2)+1))
# println("EEE = ", entanglement_entropy(mps, Int(N/2)))
# state = quantum_state_coefficients(mps_GHZ, N)
# display(state)

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking that the minimum energy from the variational ground state search agrees with the minimum energy from exact diagonalization
# and checking the entanglement entropy functions

# N = 4
# x = 4.0
# m_g_ratio = 0.5
# l_0 = 0.5 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 12
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)
# println("Minimum energy from variational ground state search: ", E_0)
# println(entanglement_entropy_inefficient(mps_ground, N))
# println(entanglement_entropy(mps_ground, N))
# matrix = mpo_to_matrix(mpo)
# display(eigvals(matrix))
# println("Minimum energy from exact diagonalization: ", minimum(eigvals(matrix)))

# ----------------------------------------------------------------------------------------------------------------------------------

# Plotting entanglement entropy vs m/g

# mg_list = LinRange(-0.7, -0.4, 10)
# ee_list = []
# N = 4
# d = 2
# D = 8
# l_0 = 0.0
# x = 1.0
# lambda = 100.0
# acc = 10^-8
# ms = 100
# for mg in mg_list
#     mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
#     _, mps, _ = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)
#     append!(ee_list, entanglement_entropy(mps, N))
# end
# plot(mg_list, ee_list)

# ----------------------------------------------------------------------------------------------------------------------------------

# Check that the penalty term enforces total charge to 0 and checking the local charge MPO

# N = 10
# x = 1.0
# m_g_ratio = 0.125
# l_0 = 0.0 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 100
# d = 2
# D = 10

# penalty_mpo = get_penalty_term_MPO(N, lambda)
# penalty_mpo_matrix = mpo_to_matrix(penalty_mpo)
# penalty_matrix_exact = get_penalty_term_matrix(N, lambda)
# display(eigvals(penalty_mpo_matrix))
# display(eigvals(penalty_matrix_exact))

# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)

# mpo_penalty = get_penalty_term_MPO(2*N, lambda)
# mps_after_penalty_mpo = act_mpo_on_mps(mpo_penalty, mps_ground)

# println(get_mpo_expectation_value(mps_ground, mpo_penalty))
# println(get_spin_half_expectation_value(2*N, mps_ground, "z")) # the total charge operator is sum_n 1/2 sigma^z_n

# charge_list = [] # stores Q_tilde_n

# for n in 1:N

#     charge_mpo = get_local_charge_MPO(N, n)
#     mps_right = act_mpo_on_mps(charge_mpo, mps_ground)
#     append!(charge_list, inner_product_MPS(mps_ground, mps_right))

# end

# display(charge_list)
# println(sum(charge_list))

# l_field = get_electric_field_configuration(N, l_0, mps_ground)
# display(l_field)

# This checks Gauss' law: l_field[i] - l_field[i-1] = Q[i] ie L_n  - L_n-1 = Q_tilde_n

# for n in 1:N
#     if n == 1
#         display(isapprox(l_field[n]-l_0, charge_list[n]))
#         display(l_field[n]-l_0)
#         display(charge_list[n])
#     else
#         display(isapprox(l_field[n]-l_field[n-1], charge_list[n]))
#         display(l_field[n]-l_field[n-1])
#         display(charge_list[n])
#     end
# end

# append!(l_field, l_0)
# display(sum(l_field))

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

# Testing reading an mps from h5 file

# N = 48
# D = 80
# mg = -0.6
# x = 10.0
# mps = h5_to_mps(2*N, D, mg, x)
# println(inner_product_MPS(mps, mps))
# N = length(mps)
# display(mps[1])
# display(mps[N])
# for i in 1:N
#     println(size(mps[i]))
# end

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing particular parameters for Schwinger model

# N = 4
# D = 8
# mg = -0.6
# x = 60.0
# ms = 100
# acc = 10^(-10)
# lambda = 100.0
# l_0 = 0.0
# d = 2
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
# E_0, mps, ns = variational_ground_state_MPS_for_saving(2*N, d, D, mpo, acc, ms)

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing the gauge_mps! function

# N = 82
# D = 80
# d = 2

# mps = initialize_MPS(2*N, d, D)
# gauge_mps!(right, mps, true, 2*N)
# display(mps)

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing calculating entropy for already saved mps

# mg = 0.125
# x = 1.0
# N = 4
# D = 8
# mps_to_entropy_save_file(mg, x, 2*N, D)

# ----------------------------------------------------------------------------------------------------------------------------------

# # https://github.com/kuehnste/mps.jl/blob/main/test/runtests.jl

# @testset "MPS Testing" begin

#     @test 

# end

# ----------------------------------------------------------------------------------------------------------------------------------

