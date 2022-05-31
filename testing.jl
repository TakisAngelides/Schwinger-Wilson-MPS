using Profile
using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
# using PyPlot
using ProfileView
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

# N = 2
# x = 1.0
# m_g_ratio = -0.125
# l_0 = 0.0 # this is theta/2pi
# lambda = 2.0
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 20
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# mpo_matrix = mpo_to_matrix(mpo)
# exact_E_0, _ = eigs(mpo_matrix, nev=1, which=:SR)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)
# println("Minimum energy from variational ground state search: ", real(E_0))
# println("Exact minimum energy:                                ", real(exact_E_0[1]))
# println(entanglement_entropy_inefficient(mps_ground, N))
# println(entanglement_entropy(mps_ground, N))
# matrix = mpo_to_matrix(mpo)
# display(eigvals(matrix))
# println("Minimum energy from exact diagonalization: ", minimum(eigvals(matrix)))

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking that the minimum energy from the variational ground state search agrees in both Q=0 enforced and not enforced

# N = 20
# x = 1.0
# m_g_ratio = -0.125
# l_0 = 0.0 # this is theta/2pi
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 20

# lambda_1 = 100.0
# mpo_1 = get_Schwinger_Wilson_MPO(N, l_0, x, lambda_1, m_g_ratio)
# E_0_1, mps_ground_1, sn_1 = variational_ground_state_MPS(2*N, d, D, mpo_1, acc, max_sweeps)

# lambda_2 = 0.0
# mpo_2 = get_Schwinger_Wilson_MPO(N, l_0, x, lambda_2, m_g_ratio)
# E_0_2, mps_ground_2, sn_2 = variational_ground_state_MPS(2*N, d, D, mpo_2, acc, max_sweeps)

# println("Minimum energy from variational ground state search with zero charge enforced: ", E_0_1)
# println("Minimum energy from variational ground state search with zero charge not enforced: ", E_0_2)

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

# Testing the chiral condensate operator mpo

# N = 4
# x = 1.0
# m_g_ratio = 0.21
# l_0 = 0.4 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 30
# d = 2
# D = 8
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)
# E_0, mps_ground, sn = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps)
# operator_mpo = get_chiral_condensate_MPO(2*N)
# mps_after_operator_mpo = act_mpo_on_mps(operator_mpo, mps_ground)
# expectation_value = inner_product_MPS(mps_ground, mps_after_operator_mpo)
# println(expectation_value)

# ----------------------------------------------------------------------------------------------------------------------------------

# Generate data for total electric field vs m/g

# mg_list = LinRange(-0.7, -0.4, 5)
# N = 40
# d = 2
# D = 60
# l_0 = 0.0
# x = 10.0
# lambda = 100.0
# acc = 10^(-8)
# ms = 100
# open("E_field_vs_mass.txt", "w") do f
#     for mg in mg_list
#         mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
#         _, mps, _ = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)
#         E = sum(get_electric_field_configuration(N, l_0, mps)) + l_0
#         EE = entanglement_entropy(mps, N)
#         write(f, "$(mg),$(E),$(EE)\n")
#     end
# end

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

# Checking time for each line in the variational algorithm

# N = 10
# x = 1.0
# m_g_ratio = 0.125
# l_0 = 0.3 # this is theta/2pi
# lambda = 100.0
# acc = 10^(-10)
# max_sweeps = 100
# d = 2
# D = 10
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, m_g_ratio)

# VSCodeServer.@profview E_0_1, mps_ground_1, sn_1 = variational_ground_state_MPS(2*N, d, D, mpo, acc, max_sweeps) # Gets a flamegraph showing the amount of time for each line

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

# Testing pseudo momentum operator mpo

# N = 10
# D = 20
# mg = 0.125
# x = 1.0
# ms = 100
# acc = 10^(-10)
# lambda = 100.0
# l_0 = 100.0
# d = 2
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
# E_0, mps, ns = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)
# momentum_mpo = get_pseudo_momentum_MPO(N, x)
# E_1, mps_1, ns_1 = variational_first_excited_MPS(2*N, d, D, mpo, acc, ms, mps, E_0)
# display(get_mpo_expectation_value(mps, momentum_mpo))
# display(get_mpo_expectation_value(mps_1, momentum_mpo))

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

# Testing getting the spin configuration of an mps and then the electric field configuration
# and also testing Gauss's law L_n - L_n-1 - Q_n = 0 for 0 external charges

# N = 4
# D = 8
# mg = 2.0
# x = 1.0
# ms = 20
# acc = 10^(-8)
# lambda = 0.0
# l_0 = 2.5
# d = 2
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
# E_0, mps, ns = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)

# electric_field_configuration = get_electric_field_configuration(l_0, mps)
# charge_configuration = get_charge_configuration(mps)

# L_list = electric_field_configuration
# Q_list = charge_configuration

# println("Checking Gauss's law")
# for n in 1:N
#     if n == 1
#         println(L_list[n]-l_0-Q_list[n])
#     else
#         println(L_list[n]-L_list[n-1]-Q_list[n])
#     end
# end
# println("Charge configuration")
# display(charge_configuration)
# println("Total charge")
# println(sum(charge_configuration))
# println("Average electric field")
# println(mean(electric_field_configuration))

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking the first order phase transition

# theta_list = LinRange(0, 2*pi, 10)

# N = 64
# D = 60
# mg = 5.0
# x = 10.0
# ms = 20
# acc = 10^(-8)
# lambda = 0.0
# d = 2
# avg_E_field_list = []

# for theta in theta_list

#     println(theta)

#     l_0 = theta/(2*pi)
#     mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
#     E_0, mps, ns = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)
#     append!(avg_E_field_list, real(mean(get_electric_field_configuration(l_0, mps))))

# end

# display(avg_E_field_list)
# plot(theta_list, avg_E_field_list)
# scatter!(theta_list, avg_E_field_list)
# savefig("avg_E_vs_theta.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

# # https://github.com/kuehnste/mps.jl/blob/main/test/runtests.jl

# @testset "MPS Testing" begin

#     @test 

# end

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking Stefan's and Takis's mpo for Schwinger Wilson spectrum

# mpo_takis =         get_Schwinger_Wilson_MPO(4, 0.2, 0.3, 10.0, 0.125)
# mpo_stefan = get_Schwinger_Wilson_MPO_Stefan(4, 0.2, 0.3, 10.0, 0.125)
# mpo_gen = get_Schwinger_Wilson_general_r_MPO(4, 0.2, 0.3, 10.0, 0.125, 1.0)

# matrix_takis = mpo_to_matrix(mpo_takis)
# matrix_stefan = mpo_to_matrix(mpo_stefan)
# matrix_gen = mpo_to_matrix(mpo_gen)

# e_takis = eigvals(matrix_takis)
# e_stefan = eigvals(matrix_stefan)
# e_gen = eigvals(matrix_gen)

# for i in 1:length(e_takis)
#     println(e_takis[i], " -- ", e_stefan[i], " -- ", e_gen[i])
# end

# display(e_takis)
# display(e_stefan)
# display(e_gen)

# ----------------------------------------------------------------------------------------------------------------------------------

# Checking if the Schwinger Wilson spectrum is the same for theta = 0 m/g = -0.125 and theta = pi m/g = 0.125

# mpo_1 = get_Schwinger_Wilson_MPO(4, 0.0, 1.0, 100.0, -0.125)
# mpo_2 = get_Schwinger_Wilson_MPO(4, 0.5, 2.0, 100.0, 0.125)

# matrix_1 = mpo_to_matrix(mpo_1)
# matrix_2 = mpo_to_matrix(mpo_2)

# display(eigvals(matrix_1))
# display(eigvals(matrix_2))
# println(eigvals(matrix_1) == eigvals(matrix_2))

# mpo_2 = get_Schwinger_Wilson_MPO_Stefan(4, 0.0, 1.0, 100.0, -0.125)
# m_2 = mpo_to_matrix(mpo_2)
# display(eigvals(m_2))

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing the KrylovKit eigsolve as well as general r = 1.0 MPO and original r = 1 MPO

# N = 10
# D = 20
# x = 1.0
# ms = 100
# acc = 10^(-8)
# lambda = 100.0
# l_0 = 0.5
# d = 2
# mg = 0.125

# mpo_gen = get_Schwinger_Wilson_general_r_MPO(N, l_0, x, lambda, mg, 1.0)
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)

# E_gen, _, _ = variational_ground_state_MPS(2*N, d, D, mpo_gen, acc, ms)
# E, _, _ = variational_ground_state_MPS(2*N, d, D, mpo, acc, ms)

# println(E_gen, " -- ", E)

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing the E.E. vs m/g for r = 1 vs r = -1

# N = 20
# D = 20
# x = 1.0
# ms = 100
# acc = 10^(-8)
# lambda = 100.0
# l_0 = 0.5
# d = 2
# m_g_list = LinRange(0.0, 7.0, 20)
# EE_pos_list = []
# EE_neg_list = []

# for mg in m_g_list
#     mpo_pos = get_Schwinger_Wilson_general_r_MPO(N, l_0, x, lambda, mg, 1.0)
#     mpo_neg = get_Schwinger_Wilson_general_r_MPO(N, l_0, x, lambda, mg, -1.0)
    
#     _, mps_pos, _ = variational_ground_state_MPS(2*N, d, D, mpo_pos, acc, ms)
#     _, mps_neg, _ = variational_ground_state_MPS(2*N, d, D, mpo_neg, acc, ms)
    
#     ee_pos = entanglement_entropy(mps_pos, N)
#     append!(EE_pos_list, ee_pos)
    
#     ee_neg = entanglement_entropy(mps_neg, N)
#     append!(EE_neg_list, ee_neg)
# end

# # PyPlot.clf()
# # figure()
# plot(m_g_list, EE_pos_list)
# savefig("rpos1.pdf")
# # gcf()
# # PyPlot.clf()
# # figure()
# plot(m_g_list,EE_neg_list)
# savefig("rneg1.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing reading the mps generated by the function variational_ground_state_algorithm

# name_of_mps = "N_2_x_1.0_D_4_l0_0.25_mg_0.125_ms_20_acc_1.0e-8_lam_10.0"
# path_to_mps = "N_2_x_1.0_D_4_l0_0.25_mg_0.125_ms_20_acc_1.0e-8_lam_10.0.h5"
# N = 4
# mps = Vector{Array{ComplexF64}}(undef, N)

# h5open(path_to_mps, "r") do fid

#     g = fid[name_of_mps]
#     for i in 1:N
#         mps[i] = read(g["$(i)"])
#     end

# end

# N = 2
# x = 1.0
# l_0 = 0.25
# lambda = 10.0
# mg = 0.125
# mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)

# println(inner_product_MPS(mps, mps))
# println(get_mpo_expectation_value(mps, mpo))

# ----------------------------------------------------------------------------------------------------------------------------------

# Testing energy(theta 1) - energy(theta 2) vs m/g for r=1 vs r=-1 as it is expected that
# r = 1 should be lower than r=-1 for any given m/g

# N = 10
# D = 10
# x = 1.0
# mg = 0.5
# lambda = 100.0
# acc = 1.0e-8
# ms = 100
# d = 2

# mpo_left = get_Schwinger_Wilson_general_r_MPO(N, 0.1, x, lambda, mg, 1.0)
# mpo_right = get_Schwinger_Wilson_general_r_MPO(N, 0.4, x, lambda, mg, 1.0)

# E_left, _, _ = variational_ground_state_MPS(2*N, d, D, mpo_left, acc, ms)
# E_right, _, _ = variational_ground_state_MPS(2*N, d, D, mpo_right, acc, ms)

# E_diff_1 = E_left - E_right

# mpo_leftm = get_Schwinger_Wilson_general_r_MPO(N, 0.1, x, lambda, mg, -1.0)
# mpo_rightm = get_Schwinger_Wilson_general_r_MPO(N, 0.4, x, lambda, mg, -1.0)

# E_leftm, _, _ = variational_ground_state_MPS(2*N, d, D, mpo_leftm, acc, ms)
# E_rightm, _, _ = variational_ground_state_MPS(2*N, d, D, mpo_rightm, acc, ms)

# E_diff_m1 = E_leftm - E_rightm

# println(E_diff_1, " -- ", E_diff_m1)

# ----------------------------------------------------------------------------------------------------------------------------------

# Plotting the Electric field vs link

# N = 64
# l_0 = 0.295
# name_of_mps = "N_64_x_10.0_D_40_l0_0.295_mg_-0.0968_ms_100_acc_1.0e-8_lam_100.0_r_1.0"
# f = h5open("N_64_x_10.0_D_40_l0_0.295_mg_-0.0968_ms_100_acc_1.0e-8_lam_100.0_r_1.0.h5", "r")
# mps_group = f[name_of_mps]
# mps = Vector{Array{ComplexF64}}(undef, 2*N)
# for i in 1:2*N
#     mps[i] = read(mps_group["$(i)"])
# end
# close(f)
# electric_field_list = real(get_electric_field_configuration(l_0, mps))
# charge_list = real(get_charge_configuration(mps))
# link_number_list = 1:length(electric_field_list)
# scatter(link_number_list, electric_field_list, label = "Electric Field")
# scatter!(link_number_list, charge_list, label = "Charge Density")
# savefig("EFnQDvsSite.pdf")

# ----------------------------------------------------------------------------------------------------------------------------------
