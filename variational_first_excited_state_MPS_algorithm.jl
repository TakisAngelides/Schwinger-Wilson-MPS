using Profile
using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
include("variational_MPS_algorithm.jl")

function initialize_left_right_states(mps_first::Vector{Array{ComplexF64}}, mps_0::Vector{Array{ComplexF64}}, N::Int64)

    p_states = Vector{Array{ComplexF64}}(undef, N+1)

    p_states[1] = ones(ComplexF64, 1, 1)
    p_states[N+1] = ones(ComplexF64, 1, 1)
    
    for i in N:-1:2

        p_states[i] = contraction(conj(mps_first[i]), (3,), mps_0[i], (3,))
        p_states[i] = contraction(p_states[i], (2,4), p_states[i+1], (1,2))

    end

    return p_states
    
end

function update_p_states!(sweep_direction::Form, p_states::Vector{Array{ComplexF64}}, M::Array{ComplexF64}, mps_0_tensor::Array{ComplexF64}, i::Int64)

    site = contraction(conj(M), (3,), mps_0_tensor, (3,))

    if sweep_direction == right # Right moving sweep from left to right
    
        p_states[i] = contraction(p_states[i-1], (1,2), site, (1,3))
    
    else # Left moving sweep from right to left

        p_states[i] = contraction(site, (2,4), p_states[i+1], (1,2))

    end
end

function get_H_projection(left, mps_0_tensor, right)

    tmp = contraction(left, (2,), mps_0_tensor, (1,)) 
    B = contraction(tmp, (2,), right, (2,)) 
    P = contraction(B, (), conj(B), ()) 
    P = permutedims(P, (2,1,3,5,4,6)) 
    dimensions = size(P)
    P = reshape(P, (dimensions[1]*dimensions[2]*dimensions[3], dimensions[4]*dimensions[5]*dimensions[6]))
    
    return P

end

function get_updated_site_first_excited(L::Array{ComplexF64}, W::Array{ComplexF64}, R::Array{ComplexF64}, left::Array{ComplexF64}, right::Array{ComplexF64}, mps_0_tensor::Array{ComplexF64}, E_0)

    H, dimensions = get_Heff(L, W, R) # H is a matrix, H_(sigma_l,a_l-1,a_l)(sigma_l_dash,a_l-1_dash,a_l_dash)
    P = get_H_projection(left, mps_0_tensor, right)
    H_eff = H - E_0*P
    E, M = eigs(H_eff, nev=1, which=:SR) # nev = 1 => it will return only 1 number of eigenvalues, SR => compute eigenvalues which have the smallest real part (ie the ground state energy and upwards depending on nev), also note M'*M = 1.0+0.0im
    M = reshape(M, (dimensions[1], dimensions[2], dimensions[3])) # M is reshaped in the form sigma_i, a_i-1, a_i
    M = permutedims(M, (2,3,1)) # M is permuted into the form a_i-1, a_i, sigma_i

    return M, E[1]

end

function variational_first_excited_MPS(N::Int64, d::Int64, D::Int64, mpo::Vector{Array{ComplexF64}}, accuracy::Float64, max_sweeps::Int64, mps_0::Vector{Array{ComplexF64}}, E_0)

    mps_first = initialize_MPS(N, d, D)
    gauge_mps!(right, mps_first, true, N) 

    states = initialize_L_R_states(mps_first, mpo, N) 
    p_states = initialize_left_right_states(mps_first, mps_0, N)
    E_initial = 10^(-5) 
    E_optimal = 0 
    sweep_number = 0 
    US = 0 

    while(true)
        
        E = 0 
        
        for i in 1:N-1 

            L = states[i]
            W = mpo[i]
            R = states[i+1]
            p_left = p_states[i]
            p_right = p_states[i+1]
            M, _ = get_updated_site_first_excited(L, W, R, p_left, p_right, mps_0[i], E_0)
            mps_first[i], _ = gauge_site(left, M)
            update_states!(right, states, mps_first[i], W, i+1) 
            update_p_states!(right, p_states, mps_first[i], mps_0[i], i+1)

        end

        for i in N:-1:2 

            L = states[i]
            W = mpo[i]
            R = states[i+1]
            p_left = p_states[i]
            p_right = p_states[i+1]
            M, E = get_updated_site_first_excited(L, W, R, p_left, p_right, mps_0[i], E_0)
            US, mps_first[i] = gauge_site(right, M) 
            update_states!(left, states, mps_first[i], W, i)
            update_p_states!(left, p_states, mps_first[i], mps_0[i], i)

        end
        
        fractional_energy_change = abs((E - E_initial)/E_initial)

        if fractional_energy_change < accuracy

            E_optimal = E

            mps_first[1] = contraction(mps_first[1], (2,), US, (1,))
            mps_first[1] = permutedims(mps_first[1], (1,3,2))

            # println("Desired accuracy reached.")

            break
        
        elseif max_sweeps < sweep_number
            
            E_optimal = E

            mps_first[1] = contraction(mps_first[1], (2,), US, (1,))
            mps_first[1] = permutedims(mps_first[1], (1,3,2))

            println("Maximum number of sweeps reached before desired accuracy.")

            break

        end

        E_initial = E
        sweep_number = sweep_number + 1
    end

    return E_optimal, mps_first, sweep_number

end