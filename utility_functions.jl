using Profile
using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
using Dates
include("MPO.jl")

@enum Form begin
    left
    right
end

function get_GHZ_mps(N::Int64)::Vector{Array{ComplexF64}}

    """
    The GHZ state is (1/sqrt(2))(|00..0> + |11..1>). This function returns the exact mps representation of the GHZ state.

    Input:

    N = number of lattice sites (Integer)

    Output:

    mps = the exact mps representation of the GHZ state (Vector of arrays of complex floats)
    """

    d = 2
    D = 2
    num = 2^(-(1/(2*N)))
    M_first = zeros(ComplexF64, 1, D, d)
    M_last = zeros(ComplexF64,D, 1, d)
    M = zeros(ComplexF64,D, D, d)
    M_first[1, 1, 1] = num
    M_first[1, 2, 2] = num
    M_last[1, 1, 1] = num
    M_last[2, 1, 2] = num
    M[1, 1, 1] = num
    M[2, 2, 2] = num

    mps = Vector{Array{ComplexF64}}(undef, N)

    for i in 1:N
        if i == 1
            mps[i] = M_first
        elseif i == N
            mps[i] = M_last
        else
            mps[i] = M
        end
    end

    return mps

end

function get_electric_field_configuration(N::Int64, l_0::Float64, mps)::Vector{Float64}
    
    """
    Gets the L_n = l_0 + sum_k=1^N Q_k for the Schwinger model

    Inputs:

    N = number of physical lattice sites

    l_0 = backgroud electric field

    mps = to calculate the electric field configuration of 

    Output:

    electric_field_list = the index of this list is the left edge of the link on which the electric field is described

    """
    
    charge_list = []
    electric_field_list = []

    for n in 1:N

        charge_mpo = get_local_charge_MPO(N, n)
        mps_right = act_mpo_on_mps(charge_mpo, mps)
        append!(charge_list, real(inner_product_MPS(mps, mps_right)))
    
    end

    for n in 1:N
        
        append!(electric_field_list, l_0 + sum(charge_list[1:n]))
    
    end

    return electric_field_list

end

function operator_to_sparse_matrix(operator, idx::Int64, N::Int64)

    """
    This function will convert an operator acting on a single site with index idx to a sparse matrix.
    It uses the kron from LinearAlgebra which is the tensor product between matrices.

    Inputs:

    N = number of lattice sites (if the physical sites are M then N = 2M)

    operator = the operator in matrix representation to be converted to a sparse matrix

    idx = the index for the site on which the operator acts on

    Output:

    result = the sparse matrix representing the operator in the full space of N sites (Matrix)

    """

    result = ones(1)
    I = [1 0; 0 1]
    
    for i in 1:N
        if i == idx
            result = kron(result, operator) # this is the kronecker tensor product between two vectors or two matrices
        else
            result = kron(result, I)
        end
    end

    return result

end

function get_Schwinger_hamiltonian_matrix(N::Int64, l_0::Float64, x::Float64, lambda::Float64, m_g_ratio::Float64)

    """
    This function returns the Schwinger model Hamiltonian with Wilson fermions in the form of a matrix.

    Inputs:

    N = number of lattice sites (if M is the number of physical lattice sites then N = 2M)

    l_0 = is the background electric field, in this case l_0 = theta/2pi

    x = 1/(a*g)^2 where a is the lattice spacing and g is the coupling constant of the theory
    
    lambda = the penalty term's lagrange multiplier

    m_g_ratio = m/g where m is the fermion mass

    Output:

    result = the Schwinger model Hamiltonian as a matrix (Matrix)
    """

    A = -2*1im*(sqrt(x)*m_g_ratio + x)
    B = -2*1im*x
    C = (N-1)*l_0^2 + lambda*N/2 + N*(N-1)/4

    Z = [1 0; 0 -1]
    PLUS = [0 1; 0 0]
    MINUS = [0 0; 1 0]

    result = zeros((2^(2*N), 2^(2*N)))

    for n in 1:N

        result += A.*(operator_to_sparse_matrix(PLUS, 2*n-1, 2*N)*operator_to_sparse_matrix(MINUS, 2*n, 2*N))
        result += - A.*(operator_to_sparse_matrix(MINUS, 2*n-1, 2*N)*operator_to_sparse_matrix(PLUS, 2*n, 2*N))

    end

    for n in 1:N-1
    
        result += B.*(operator_to_sparse_matrix(MINUS, 2*n-1, 2*N)*operator_to_sparse_matrix(Z, 2*n, 2*N)operator_to_sparse_matrix(Z, 2*n+1, 2*N)operator_to_sparse_matrix(PLUS, 2*n+2, 2*N))
        result += -B.*(operator_to_sparse_matrix(PLUS, 2*n-1, 2*N)*operator_to_sparse_matrix(Z, 2*n, 2*N)operator_to_sparse_matrix(Z, 2*n+1, 2*N)operator_to_sparse_matrix(MINUS, 2*n+2, 2*N))
    
    end

    for k in 1:N-1

        result += (l_0*(N-k)).*(operator_to_sparse_matrix(Z,2*k-1,2*N) + operator_to_sparse_matrix(Z,2*k,2*N))

    end

    for k in 1:N-1
        for q in k+1:N
        
            result += (0.5*(N-q+lambda)).*(operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*q-1,2*N)+operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*q,2*N)+operator_to_sparse_matrix(Z,2*k,2*N)*operator_to_sparse_matrix(Z,2*q-1,2*N)+operator_to_sparse_matrix(Z,2*k,2*N)*operator_to_sparse_matrix(Z,2*q,2*N))

        end
    end

    for k in 1:N
    
        result = result + (0.5*(N-k+lambda)).*(operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*k,2*N))
    
    end

    result += C.*diagm(ones(2^(2*N)))

    return result

end

function get_penalty_term_matrix(N::Int64, lambda::Float64)

    """
    This will return the penalty term of the W' dimensionless Hamiltonian in the form of a matrix.
    
    Inputs:

    N = number of lattice sites (if M is the number of physical lattice sites then N = 2M)

    lambda = the lagrange multiplier to the penalty term in W'

    Output:

    result = the penalty term in the form of a matrix (Matrix)
    """

    Z = [1 0; 0 -1]

    result = zeros((2^(2*N), 2^(2*N)))

    for k in 1:N-1
        for q in k+1:N
        
            result += (0.5*(lambda)).*(operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*q-1,2*N)+operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*q,2*N)+operator_to_sparse_matrix(Z,2*k,2*N)*operator_to_sparse_matrix(Z,2*q-1,2*N)+operator_to_sparse_matrix(Z,2*k,2*N)*operator_to_sparse_matrix(Z,2*q,2*N))

        end
    end

    for k in 1:N
    
        result = result + (0.5*(lambda)).*(operator_to_sparse_matrix(Z,2*k-1,2*N)*operator_to_sparse_matrix(Z,2*k,2*N))
    
    end

    result += (lambda*N/2).*diagm(ones(2^(2*N)))

    return result

end

function mpo_to_matrix(mpo::Vector{Array{ComplexF64}})
    """
    Converts an MPO to a matrix by contracting the bond links and reshaping into (sigma_1...sigma_N),(sigma'_1...sigma'_N)
    
    Input:

    mpo = the mpo to be converted to a matrix (Vector of arrays of complex numbers)

    Output:

    result = the matrix with indices (sigma_1...sigma_N),(sigma'_1...sigma'_N) 
    """

    N = length(mpo)

    idx = 2

    result = contraction(mpo[1], (idx,), mpo[2], (1,))

    for i in 3:N
        
        idx = idx + 2
        result = contraction(result, (idx,), mpo[i], (1,))

    end
        
    result = contraction(result, (length(size(result))-2,), ones(1), (1,))
    result = contraction(result, (1,), ones(1), (1,))

    odd = 1:2:2*N-1
    even = 2:2:2*N

    result = permutedims(result, (odd..., even...))

    dims = size(result)
    total_dim = dims[1:N]
    total_dim_dash = dims[N+1:2*N]

    result = reshape(result, (prod(total_dim), prod(total_dim_dash)))

    return result
end

function act_mpo_on_mps(mpo::Vector{Array{ComplexF64}}, mps::Vector{Array{ComplexF64}})::Vector{Array{ComplexF64}}

    """
    Act with an mpo on an mps to produce a new mps with increased bond dimension.

    Inputs:

    mpo = the mpo to act on the mps (Vector of Arrays)

    mps = the mps that the mpo will act on (Vector of Arrays)

    Output:

    result = the new mps with increased bond dimension resulting from acting with the mpo input on the mps input
    """
    
    N = length(mps)

    result = Vector{Array{ComplexF64}}(undef, N)

    for i in 1:N
    
        tmp = contraction(mpo[i], (4,), mps[i], (3,)) # Does contraction of sigma'_i: W_(b_i-1 b_i sigma_i sigma'_i) M_(a_i-1 a_i sigma'_i) = T_(b_i-1 b_i sigma_i a_i-1 a_i)
        tmp = permutedims(tmp, (4, 1, 5, 2, 3)) # T_(a_i-1 b_i-1 a_i b_i sigma_i)
        idx_dims = size(tmp) # returns a tuple of the dimensions of the indices of the tensor T_(a_i-1 b_i-1 a_i b_i sigma_i)
        result[i] = reshape(tmp, (idx_dims[1]*idx_dims[2], idx_dims[3]*idx_dims[4], idx_dims[5])) # merges the bond indices of i-1 together and the bond indices of i together by reshaping the tensor into having indices of higher dimension 

    end

    return result

end

function get_spin_half_expectation_value(N::Int64, mps::Vector{Array{ComplexF64}}, measure_axis::String)::ComplexF64

    """
    Computes the expectation value of the operator

    O^j = sum_i (S^j_i) = sum_i (sigma^j_i/2)

    which gives the total magnetisation along a specific axis j which can be x, y or z and i runs over all sites from 1 to N.

    Inputs:

    N = number of lattice sites (Integer)

    mps = the mps which we will use to calculate <mps|operator|mps> (Vector of 3-tensors)

    measure_axis = which axis to measure the spin (String) takes values "x", "y" or "z"

    Outputs:

    result = total magnetisation of spin 1/2 along a given x, y or z axis (ComplexF64)
    """

    # If we want to measure spin in the x direction we get the MPO operator = sum_i sigma^x_i/2 or y or z equivalently

    if measure_axis == "x"
        mpo = get_spin_half_MPO(N, "x")
    elseif measure_axis == "y"
        mpo = get_spin_half_MPO(N, "y")
    else
        mpo = get_spin_half_MPO(N, "z")
    end

    # Contracts the triple of <mps|mpo|mps> at site 1, then contracts this triple with a dummy 1x1x1 tensor of value 1
    # which will get rid of the trivial indices of the first triple at site 1. The trivial indices are the ones labelled 1,
    # see for example Schollwock equation (192) first bracket.

    triple_1 = contraction(conj!(deepcopy(mps[1])), (3,), mpo[1], (3,))
    triple_1 = contraction(triple_1, (5,), mps[1], (3,))
    dummy_tensor = ones(ComplexF64, 1,1,1)
    result = contraction(dummy_tensor, (1,2,3), triple_1, (1,3,5))

    # Now we compute the triple <mps|mpo|mps> at site i and contract it with the triple at site i-1 which we named result before

    for i in 2:N
    
        triple_i = contraction(conj!(deepcopy(mps[i])), (3,), mpo[i], (3,))
        triple_i = contraction(triple_i, (5,), mps[i], (3,))
        result = contraction(result, (1,2,3), triple_i, (1,3,5))

    end

    return result[1,1,1] # expectation value of total magnetisation with respect to a give x,y or z axis which was a 1x1x1 tensor hence the [1,1,1] index to get a ComplexF64

end

function get_mpo_expectation_value(mps::Vector{Array{ComplexF64}}, mpo::Vector{Array})::ComplexF64

    """
    Computes the expectation value of the operator represented by the mpo input

    Inputs:

    mps = the mps which we will use to calculate <mps|operator|mps> (Vector of 3-tensors)

    mpo = the mpo representing the operator to evaluate the expectation value on

    Outputs:

    result = expectation value of mpo with respect to mps
    """

    N = length(mps)

    @assert(length(mpo) == N, "The length of the mps is $(N) and the length of the mpo is $(length(mpo)), hence cannot take the expectation value of the mpo.")

    result = ones(Float64, 1, 1, 1)

    for i in 1:N
        result = contraction(result, (1,), conj(mps[i]), (1,))
        result = contraction(result, (1, 4), mpo[i], (1, 3))
        result = contraction(result, (1, 4), mps[i], (1, 3))
    end

    return result[1,1,1]

end

function inner_product_MPS(mps_1::Vector{Array{ComplexF64}}, mps_2::Vector{Array{ComplexF64}})::ComplexF64

    """
    Computes the inner product of two MPS as <mps_1|mps_2>

    Note 1: See Schollwock equation (95)

    Note 2: See my personal website for notes on this function and how the contractions are done in a specific order for efficiency.

    Inputs:

    mps_1 = The bra MPS state of the inner product (Vector)

    mps_2 = The ket MPS state of the inner product (Vector)

    Output:

    result = The complex value of the inner product (Complex)
    """

    # Assert that the number of sites in each MPS are equal

    N = length(mps_1)

    @assert(N == length(mps_2), "The two MPS inputs do not have the same number of sites.")
 
    # conj takes the input to its complex conjugate not complex transpose. Also note the first contraction which happens 
    # at the very left contracts the first index of the bra matrix with the first index of the ket matrix but these indices are 
    # trivial and set to 1. It also contracts the physical indices of the two aforementioned matrices which gives a new non-trivial result.

    result = contraction(conj(mps_1[1]), (1, 3), mps_2[1], (1, 3))
    
    for i in 2:N

        result = contraction(result, (2,), mps_2[i], (1,))

        result = contraction(conj(mps_1[i]), (1, 3), result, (1, 3))
    
    end
    
    return result[1] # results ends up being a 1x1 matrix that is why we index it with [1] to get its value
        
end

function quantum_state_coefficients(mps::Vector{Array{ComplexF64}}, N::Int64)::Array{ComplexF64}

    """
    If we write a quantum state as psi = psi_sigma1,sigma2...|sigma1>|sigma2>... then this function returns the tensor
    psi_sigma1,sigma2 of the psi represented by the input MPS.

    Inputs:

    mps = the mps that represents the quantum state for which we want the coefficients (Vector with elements being 3-tensors ie 3-arrays)

    N = number of lattice sites (Integer)

    Outputs:

    result = coefficients of quantum state namely the psi_sigma1,sigma2,... coefficients (Array of complex floats 64)

    """

    result = contraction(mps[1], (2,), mps[2], (1,))
    for i in 2:N-1
        result = contraction(result, (i+1,), mps[i+1], (1,))
    end

    result = contraction(ones(ComplexF64, 1), (1,), result, (1,))
    result = contraction(ones(ComplexF64, 1), (1,), result, (N,))

    return result

end

function gauge_site_center_orthogonal(form::Form, M_initial::Array)::Tuple{Array{ComplexF64}, Vector{Float64}, Matrix{ComplexF64}}

    """
    Gauges a site into left or right canonical form

    Note 1: See Schollwock equations (136), (137) at link: https://arxiv.org/pdf/1008.3477.pdf

    Inputs: 

    form = left or right depending on whether we want the site in left or right canonical form (of enumarative type Form)

    M_initial = 3-array to gauge representing the 3-tensor on a given site (Array)

    Output:

    If input form is left: A, S, Vt # A_(a_i-1)(s_i)(sigma_i), S_(s_i)(s_i), Vt_(s_i)(a_i)

    If input form is right: U, S, B # U(a_i-1)(s_i-1), S_(s_i-1)(s_i-1), B_(s_i-1)(a_i)(sigma_i)
    """
    
    # Julia is call by reference for arrays which are mutable so manipulations on M_initial in this function will reflect on the original unless we remove that reference with eg M = permutedims(M_initial, (1,2,3))

    if form == right # See Schollwock equation (137) for right canonical form (link: https://arxiv.org/pdf/1008.3477.pdf)

        D_left, D_right, d = size(M_initial) # Dimensions of indices of site represented by M_initial to be SVD decomposed
        # The next line is enough to remove the reference on M_initial so that it does not mutate the original M_initial and just uses its value, hence the gauge_site function does not mutate M_initial at all
        M = permutedims(M_initial, (1,3,2)) # Assumes initial index was left right physical and now M_(a_i-1)(sigma_i)(a_i)
        M = reshape(M, (D_left, d*D_right)) # Merging indices: Prepare as 2 index tensor to give to SVD, M_(a_i-1)(sigma_i)(a_i) -> M_(a_i-1)(sigma_i a_i)
        F = svd(M) # One can recover M by M = U*Diagonal(S)*Vt 
        U = F.U # U_(a_i-1)(s_i-1)
        S = F.S # S_(s_i-1)(s_i-1) although S here is just a vector storing the diagonal elements
        
        # @test length(S) == min(D_left, d*D_right) # property of SVD, note S is returned as a vector not a diagonal matrix

        # Note for complex M_initial, the following should be named Vd for V_dagger rather than Vt for V_transpose but we keep it Vt
        Vt = F.Vt # Vt_(s_i-1)(sigma_i a_i)
        Vt = reshape(Vt, (length(S), d, D_right)) # Unmerging indices: Vt_(s_i-1)(sigma_i a_i) -> Vt_(s_i-1)(sigma_i)(a_i)
        B = permutedims(Vt, (1,3,2)) # Vt_(s_i-1)(sigma_i)(a_i) -> B_(s_i-1)(a_i)(sigma_i)

        # @test isapprox(contraction(B, (2,3), conj(B), (2,3)), I) # right canonical form property

        return U, S, B # US_(a_i-1)(s_i-1), B_(s_i-1)(a_i)(sigma_i)

    else # See Schollwock equation (136) for left canonical form

        D_left, D_right, d = size(M_initial)
        M = permutedims(M_initial, (3, 1, 2)) # M_(a_i-1)(a_i)(sigma_i) -> M_(sigma_i)(a_i-1)(a_i)
        M = reshape(M, (d*D_left, D_right)) # M_(sigma_i)(a_i-1)(a_i) -> M_(sigma_i a_i-1)(a_i)
        F = svd(M)
        U = F.U # U_(sigma_i a_i-1)(s_i)
        S = F.S # S_(s_i)(s_i) although stored as vector here

        # @test length(S) == min(d*D_left, D_right) # property of SVD, note S is returned as a vector not a diagonal matrix

        Vt = F.Vt # Vt_(s_i)(a_i)
        U = reshape(U, (d, D_left, length(S))) # U_(sigma_i)(a_i-1)(s_i)
        A = permutedims(U, (2, 3, 1)) # A_(a_i-1)(s_i)(sigma_i)

        # @test isapprox(contraction(conj(A), (1,3), A, (1,3)), I) # left canonical form property

        return A, S, Vt # A_(a_i-1)(s_i)(sigma_i), SVt_(s_i)(a_i)

    end

end

function center_orthogonalize!(mps::Vector{Array{ComplexF64}}, idx::Int64)::Vector{Float64}

    """
    This function puts/mutates an mps in center orthogonal form about the site with index idx so that everything to the left of this idx
    including the idx itself is in left canonical form and everything to the right is in right canonical form.

    We start from site 1 and put everything up to idx in left canonical form. When we gauge the site at index idx we are left with 
    the matrix SVt. Then we start from site N (the last site of the mps) and put everything in right canonical form up to and including
    the site at index idx+1. When we gauge the site at index idx+1 we are left with the matrix US. We then compute SVt*US to get a 
    matrix which we call M_tilde. We perform SVD on this M_tilde to get U, S, Vt. U is multiplied on the site at index idx and 
    Vt is multiplied on the site at index idx+1. Note that these latter multiplications preserve the left and right canonical forms
    at each site. We return S which stores the singular values of this bi-partition of the mps. The end result looks schematically like
    L-L-S-R-R for N = 4 and idx = 2 where the initial mps would have looked like M-M-M-M.

    Note 1: 
    
    This function preserves the norm of the mps if one computes the norm using the function norm_center_orthogonal. If one uses the
    function inner_product_MPS then that would give norm 1, ie a normalised mps, as expected since everything is in left/right canonical form.

    Note 2:

    This function returns the singular values in the form of a vector of floats. If one wants to convert that to a matrix one can call Diagonal(S).

    Inputs:

    mps = the mps to gauge into a center orthogonal form at index idx (Vector of arrays of complex floats)

    idx = the index of the site for which everything to the right is in right canonical form and everything to the left including the site at index idx is in left canonical form (Integer)

    Output:

    S = the singular values which reside at the link between sites at index idx and idx+1 (Vector of floats)
    """
    
    # Put everything from site 1 up to and including site at index idx in left canonical form and save the last SVt matrix left over from
    # left gauging the tensor at site with index idx.

    M_tilde = mps[1]
    SVt = undef
    for i in 1:idx
    
        mps[i], SVt = gauge_site(left, M_tilde)
        M_tilde = contraction(SVt, (2,), mps[i+1], (1,)) # multiply SVt left over matrix on the right neighbour
    
    end

    # Put everything from site N down to and including site at index idx+1 in right canonical form and save the last US matrix left over from
    # right gauging the tensor at site with index idx+1.

    N = length(mps)
    M_tilde = mps[N]
    US = undef

    for i in N:-1:idx+1
        
        US, mps[i] = gauge_site(right, M_tilde)
        M_tilde = contraction(mps[i-1], (2,), US, (1,)) # mutliply US left over matrix on the left neighbour
        M_tilde = permutedims(M_tilde, (1,3,2)) # so as to keep the physical index last 
    
    end

    # Calculate M_tilde = SVt*US and perform SVD on it

    M_tilde = contraction(SVt, (2,), US, (1,))
    
    F = svd(M_tilde)
    U = F.U
    S = F.S
    B = F.Vt

    # Multiply U on idx and Vt (which we also call B) on idx+1 

    mps[idx] = contraction(mps[idx], (2,), U, (1,)) 
    mps[idx+1] = contraction(B, (2,), mps[idx+1], (1,))
    mps[idx] = permutedims(mps[idx], (1,3,2))
    
    return S

end

function norm_center_orthogonal(mps::Vector{Array{ComplexF64}}, S::Vector{Float64}, idx::Int64)::ComplexF64

    """
    This function complements the function center_orthogonalize! and returns the norm of an mps put in center orthogonal form.
    The only difference from the function inner_product_MPS is that this function multiplies the matrix of singular values
    residing at the link between index idx and idx+1 onto the tensor residing at the site with index idx.

    Inputs:

    mps = the mps to return the norm of (Vector of arrays of complex floats)

    S = the singular values residing at the link between idx and idx+1 sites (Vector of floats)

    idx = the index to the last site which is in left canonical form (Integer)

    Output:

    result[1] = the norm of the input mps (complex float)
    """

    N = length(mps)

    M_tilde = contraction(mps[idx], (2,), Diagonal(S), (1,)) # multiply the singular values on the tensor at index idx
    mps[idx] = permutedims(M_tilde, (1, 3, 2)) # so as to keep the physical index last

    # Perform standard norm calculation - see the function inner_product_MPS
    
    result = contraction(conj(mps[1]), (1, 3), mps[1], (1, 3))
    
    for i in 2:N

        result = contraction(result, (2,), mps[i], (1,))
        result = contraction(conj(mps[i]), (1, 3), result, (1, 3))
    
    end
    
    return result[1] 

end

function entanglement_entropy_inefficient(mps::Vector{Array{ComplexF64}}, idx::Int64)::Float64

    """
    Gets the entanglement entropy for the bi-partition of an mps given as input where the idx input specifies the first site
    belonging to the right partition. The right partition consists of all sites from idx until the last site and everything else
    belongs to the left partition. 

    The way this function works is by contracting all the bond links for |psi> and <psi| which together form the density matrix. 
    Then it contracts the physical indices of the right partition to get the reduced density matrix for the left partition.
    Then we calculate the eigenvalues of the reduced density matrix after we reshape it to a matrix and compute the entanglement
    entropy using the standard formula.

    Note:

    This function is highly memory inefficient since the result of contracting the bond indices leaves us with a tensor which has
    N = lenght(mps) physical indices resultin in d^N components where d is the number of physical degrees of freedom. Hence this 
    function should only be used for cross-checking other more memory efficient functions.

    Inputs:

    mps = the state for which we find the density matrix, get the reduced matrix for the left partition and calcuate its entanglement entropy (Vector of arrays of complex floats)

    idx = the index to the first site of the right partition (Integer)

    Outputs:

    S_left = the entanglement entropy of the left partition (float)
    """

    N = length(mps)
    d = size(mps[1])[3]

    @assert(idx >= 2 && idx <= N, "The index that defines the start of the right partition for the entanglement entropy
    should satisfy the inequality 2<=idx<=N, with the first equality implying that the left partition consists only of the first left
    most site and the second equality implying that the left partition consists of only sites except the last right
    most site. The input was idx = $(idx), which does not satisfy the inequality.")

    gauge_mps!(right, mps, true, N) # we normalise the mps just in case it is not given in normalised otherwise this function's output will be wrong

    tmp = ones(ComplexF64, 1)
    
    # Contract all the bond links for |psi>

    psi_bond_contracted = contraction(mps[N], (2,), tmp, (1,)) # result has 1 bond index and 1 physical index
    for i in N-1:-1:1
        psi_bond_contracted = contraction(mps[i], (2,), psi_bond_contracted, (1,)) # each time we add 1 physical index and contract two bond indices (this is done N-1 times)
    end
    psi_bond_contracted = contraction(tmp, (1,), psi_bond_contracted, (1,)) # here we just remove a trivial index, result has N bond indices and N physical indices
    
    # Contract all the bond links for <psi|

    psi_bar_bond_contracted = contraction(conj(mps[N]), (2,), tmp, (1,)) 
    for i in N-1:-1:1
        psi_bar_bond_contracted = contraction(conj(mps[i]), (2,), psi_bar_bond_contracted, (1,))
    end
    psi_bar_bond_contracted = contraction(tmp, (1,), psi_bar_bond_contracted, (1,))

    # Get the indices corresponding to the sites that belong to the right partition which will be contracted along the physical indices
    # between |psi> and <psi|

    nums = [i for i in idx:N] # indices of the right partition
    nums = tuple(nums...)
    rho_left_reduced = contraction(psi_bond_contracted, nums, psi_bar_bond_contracted, nums) # tracing out the right partition
    rho_left_reduced = reshape(rho_left_reduced, (d^(idx-1), d^(idx-1))) # reduced density matrix for the left partition after tracing out the right partition

    evals = eigvals(rho_left_reduced) # eigenvalues of the reduced density matrix
    evals = abs.(real(evals))
    S_left = 0 # will store the entanglement entropy coming from the reduced density matrix of the left partition
    for eval in evals
        if eval > 10^(-12)
            S_left += -eval*log2(eval)
        end
    end

    # For the right partition reduced density matrix - note S_right should be the same as S_left

    # nums = [i for i in 1:idx-1]
    # nums = tuple(nums...)
    # rho_right_reduced = contraction(psi_bond_contracted, nums, psi_bar_bond_contracted, nums)
    # rho_right_reduced = reshape(rho_right_reduced, (d^(N-idx+1), d^(N-idx+1)))

    # evals = eigvals(rho_right_reduced)
    # evals = abs.(real(evals))
    # S_right = 0
    # for eval in evals
    #     if eval > 10^(-12)
    #         S_right += -eval*log2(eval)
    #     end
    # end

    return S_left

end

function entanglement_entropy_old(mps::Vector{Array{ComplexF64}}, idx::Int64)::Float64

    """
    Calculates the entanglement entropy of the mps bi partitioned at site idx where site idx is the last site of the left partition.

    This function works by first normalising the input mps in case it is not (the function's output would be wrong otherwise) and 
    then puts everything to the left of idx in left canonical form. The site idx is the last site we gauge in the left canonical form
    which leaves us with matrices S, Vt left over. The S are the singular values we use to compute the entropy. This function works
    the same as if we are putting the mps in center orthogonal form and using the resulting S matrix at the bond separating left and 
    right orthogonal matrices.

    Note 1:

    This function works because we first gauge the whole input mps in right canonical form. The last SVD we do results in left over
    matrices S and B from the site at index idx. The S is what we use for the entanglement entropy and the B can be multiplied on 
    the right neighbour at idx+1 which is already in right canonical form and this multiplication will not change the right canonical
    form property of that site. So schematically the process for N = 4 and idx = 2 is MMMM -> BBBB -> UUS(BB)B. Remember that U are 
    left canonical and B are right canonical.

    Note 2:

    This is essentially the same as putting the mps in center orthogonal form using the function center_orthogonalize! and then
    using the resulting S matrix residing in the link of the orthogonal center.

    Inputs:

    mps = the mps state for which to calculate the entanglement entropy (Vector of arrays of complex floats)
    
    idx = the index to the last site of the left part of the bi partition of the mps (Integer)

    Output:

    ee = entanglement entropy for the bi-partition of the mps at site idx (float)
    """
    
    gauge_mps!(right, mps, true, length(mps)) # normalise the input mps in case it was not - function's output would be wrong if we dont put the mps in right canonical form here

    # Put all sites up to and including site at idx in left canonical form and save the left over S matrix from left gauging site idx
    
    M_tilde = mps[1]
    S_left_over = undef

    for i in 1:idx
    
        _, S, Vt = gauge_site_center_orthogonal(left, M_tilde)
        SVt = Diagonal(S)*Vt 
        M_tilde = contraction(SVt, (2,), mps[i+1], (1,))
        if i == idx
            S_left_over = S
        end

    end

    # Use the S matrix which is actually represented as a vector of floats to calculate the entanglement entropy with the standard formula

    evals = S_left_over.^2
    evals = abs.(real(evals))
    ee = 0
    for eval in evals
        if eval > 10^(-12)
            ee += -eval*log2(eval)
        end
    end

    return ee

end

function entanglement_entropy(mps::Vector{Array{ComplexF64}}, idx::Int64)::Float64

    """
    Gets the entanglement entropy for the bi-partition of the normalised mps at index idx where idx is the last site of the left partition.

    The function puts everything in left canonical form from site 1 up to and including site idx and saves the left over SVt matrix.
    Then puts everything in right canonical form from site N down to and cinluding site idx+1 and saves the left over US matrix.
    We then computer M_tilde = SVtUS and perform SVD on M_tilde to get the S matrix which we then use for the standard formula of
    entanglement entropy.

    Note:

    This function expects the input mps to be normalised otherwise the output will be wrong.

    Inputs:

    mps = the normalised mps state for which to calculate the entanglement entropy (Vector of arrays of complex floats)
    
    idx = the index to the last site of the left part of the bi partition of the mps (Integer)

    Output:

    ee = entanglement entropy for the bi-partition of the mps at site idx (float)
    """

    # Put sites 1 to idx in left canonical form and save the left over SVt matrix.

    M_tilde = mps[1]
    SVt = undef
    for i in 1:idx
    
        _, SVt = gauge_site(left, M_tilde)
        M_tilde = contraction(SVt, (2,), mps[i+1], (1,))
    
    end

    # Put sites idx+1 to N in right canonical form and save the left over US matrix.

    N = length(mps)
    M_tilde = mps[N]
    US = undef

    for i in N:-1:idx+1
        
        US, _ = gauge_site(right, M_tilde)
        M_tilde = contraction(mps[i-1], (2,), US, (1,))
        M_tilde = permutedims(M_tilde, (1,3,2))
    
    end

    # Compute M_tilde = SVt*US and perform SVD on M_tilde to use the resulting S matrix to calculate the entanglement entropy

    M_tilde = contraction(SVt, (2,), US, (1,))
    
    S = svdvals(M_tilde)

    evals = S.^2
    evals = abs.(real(evals))
    ee = 0
    for eval in evals
        if eval > 10^(-12)
            ee += -eval*log2(eval)
        end
    end

    return ee

end

function generate_Schwinger_data(mg, x, N, D, accuracy, lambda, l_0, max_sweep_number, D_previous)

    """
    Generates data for the Schwinger model with Wilson fermions in order to extrapolate to the continuum. It saves to an h5 file
    the ground state found for the given input parameters and it also writes to a text file the ground state energy and the
    expectation value for the penalty term just to ensure that is indeed 0 and does not affect the value of the ground state energy.

    Inputs:

    mg = m/g

    x = 1/(a*g)^2

    N = number of physical lattice sites

    D = bond dimension of to be initialized for the variational ground state search algorithm

    accuracy = the accuracy for when to stop the variational ground state search algorithm compared to the change in energy found by the algorithm

    lambda = lagrange multiplier to the penalty term

    l_0 = background electric field, this is theta/2pi

    max_sweep_number = the maximum number of sweeps to do in the variational ground state search algorithm

    D_previous = if this is 0 then we do not use a previous ground state solution with smaller D for the next search with greater D otherwise we do so

    Output:

    No output other than writing to h5 and txt files
    """

    penalty_term_MPO = get_penalty_term_MPO(2*N, lambda)

    open("observables_$(N)_$(D)_$(mg)_$(x).txt", "w") do file
        h5open("mps_$(N)_$(D)_$(mg)_$(x).h5", "w") do fid

            if D_previous != 0

                f = h5open("mps_$(N)_$(D_previous)_$(mg)_$(x).h5", "r")

                mps_group = f["$(lambda)_$(l_0)_$(mg)_$(x)_$(N)_$(D_previous)"]

                mps_previous = Vector{Array{ComplexF64}}(undef, 2*N)
                
                for i in 1:2*N
                
                    mps_previous[i] = read(mps_group["mps_$(i)"])

                end

                mps_after = Vector{Array{ComplexF64}}(undef, 2*N)

                for i in 1:2*N
                    
                    dims = size(mps_previous[i])
                    D_left_previous = dims[1]
                    D_right_previous = dims[2]

                    if i == 1

                        mps_after[i] = zeros(ComplexF64, 1, D, 2)
                        for j in 1:2
                            mps_after[i][1:1, 1:D_right_previous, j] = mps_previous[i][:, :, j]
                        end

                    elseif i == 2*N

                        mps_after[i] = zeros(ComplexF64, D, 1, 2)
                        for j in 1:2
                            mps_after[i][1:D_left_previous, 1:1, j] = mps_previous[i][:,:,j]
                        end

                    else

                        mps_after[i] = zeros(ComplexF64, D, D, 2)
                        for j in 1:2
                            mps_after[i][1:D_left_previous, 1:D_right_previous, j] = mps_previous[i][:,:,j]
                        end
                    end
                end
            else
                mps_after = 0
            end
            
            t1 = Dates.now()
            println("Now calculating: lambda = $(lambda), l_0 = $(l_0), m_over_g = $(mg), x = $(x), N = $(N), D = $(D), and the time is $(t1)\n")
            
            mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
            if D_previous != 0
                E_0, mps, sweeps = variational_ground_state_MPS_from_previous(2*N, 2, D, mpo, accuracy, max_sweep_number, D_previous, mps_after)
            else
                E_0, mps, sweeps = variational_ground_state_MPS(2*N, 2, D, mpo, accuracy, max_sweep_number)
            end
            
            if !(haskey(fid, "$(lambda)_$(l_0)_$(mg)_$(x)_$(N)_$(D)"))
                create_group(fid, "$(lambda)_$(l_0)_$(mg)_$(x)_$(N)_$(D)")
            end
            
            g = fid["$(lambda)_$(l_0)_$(mg)_$(x)_$(N)_$(D)"]
            
            for i in 1:length(mps)
                
                g["mps_$(i)"] = mps[i]
                
            end
                
            penalty_term_expectation_value = get_mpo_expectation_value(mps, penalty_term_MPO)
            
            write(file, "$(lambda),$(l_0),$(mg),$(x),$(N),$(D),$(sweeps),$(penalty_term_expectation_value),$(E_0)\n")
            
            t2 = Dates.now()
            println("Have just finished calculating: lambda = $(lambda), l_0 = $(l_0), m_over_g = $(mg), x = $(x), N = $(N), D = $(D), and the time is $(t2)\n")
            
        end
    end
end

function generate_entropy_data_Ising()

    N = 100
    D = 40
    accuracy = 10^(-10)
    max_sweep_number = 100
    J = -1.0
    g_z = -0.1
    l, u, n = 0.0, 3.0, 20
    g_x_list = LinRange(l, u, n)

    open("entropy_mass_data_Ising_$(l)_$(u)_$(n).txt", "w") do file

        for g_x in g_x_list

            mpo = get_Ising_MPO(N, J, g_x, g_z)
            E_0, mps, sweeps = variational_ground_state_MPS(N, 2, D, mpo, accuracy, max_sweep_number)
            ee = entanglement_entropy(mps, Int(N/2))   
            write(file, "$(g_x),$(ee)\n")
            
        end

    end 

end

function generate_entropy_data(mg, x, N, D, accuracy, lambda, l_0, max_sweep_number)

    mpo = get_Schwinger_Wilson_MPO(N, l_0, x, lambda, mg)
    from_saved_mps = false
    _, _, _ = variational_ground_state_MPS_for_saving(2*N, 2, D, mpo, accuracy, max_sweep_number, from_saved_mps)

end

function h5_to_mps(N::Int64, D::Int64, mg::Float64, x::Float64)::Vector{Array{ComplexF64}}

    """
    Input:

    N = the number of spin sites which is double the physical sites (Int)

    D = bond dimension (Int)

    mg = mass over coupling constant (Float)

    x = 1/(a^2 * g^2) (Float)

    Output:

    mps = the mps saved in the h5 file with the name "mps_$(N)_$(D)_$(mg)_$(x).h5"
    """

    name_of_file = "mps_$(N)_$(D)_$(mg)_$(x).h5"

    f = h5open(name_of_file, "r")

    lambda = 100.0
    l_0 = 0.0
    
    mps_group = f["$(lambda)_$(l_0)_$(mg)_$(x)_$(N)_$(D)"]

    mps = Vector{Array{ComplexF64}}(undef, N)
    
    for i in 1:N
    
        mps[i] = read(mps_group["mps_$(i)"])

    end

    close(f)
    
    return mps

end

function mps_to_entropy_save_file(mg::Float64, x::Float64, N::Int64, D::Int64)

    mps = h5_to_mps(N, D, mg, x)
    half = Int(N/2)
    ee = entanglement_entropy(mps, half)
    open("entropy_mass_data_Schwinger_$(mg)_$(x)_$(N)_$(D).txt", "w") do file
        write(file, "$(mg),$(ee)\n")
    end 

end