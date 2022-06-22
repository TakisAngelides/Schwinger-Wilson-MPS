using Profile
using LinearAlgebra
using Arpack
using BenchmarkTools
using Plots
using LaTeXStrings
using Test
using HDF5
using Dates

@enum Form begin
    left
    right
end

function get_Ising_MPO(N::Int64, J::Float64, g_x::Float64, g_z::Float64)::Vector{Array{ComplexF64}}

    """
    Creates the 1D transverse and longitudinal applied magnetic field Ising model with open boundary conditions with the Hamiltonian 
    operator 
    
    H = sum_over_nearest_neighbours(J * Z_i * Z_j) + sum_over_all_sites((g_x * X_i) + (g_z * Z_i)). 
    
    It stores the MPO as a vector of N elements, 1 element for each lattice site, and at site or in other words each element is a 
    4-tensor or in memory a 4-array with the indices stored as shown below.

          sigma_i                                              3
             |                                                 |
    alpha -- Wi -- beta , which is stored in the order: 1 -- mpo[i] -- 2
             |                                                 |
        sigma_i_dash                                           4      

    Note 1: The indices are stored in the order alpha, beta, sigma_i, sigma_i_dash. The first two are the left and right bond indices
    and the last two are the physical indices that would connect to the physical indices of the bra (above MPO) and ket (below MPO)
    respectively.

    Note 2: See notes for this function, including the derivation of the W1, WN and Wi, in Constructing Hamiltonian MPO note on my website.

    Note 3: See equation (8) at https://arxiv.org/pdf/1012.0653.pdf.

    Inputs:

    N = Number of physical sites on the lattice (Integer)
    J, g_x, g_z = coupling constants in the Hamiltonian (Floats)

    Output:

    mpo = vector that stores 4-arrays representing the MPO
    """

    D = 3 # Bond dimension for the 1D transverse and longitudinal field Ising model is 3

    d = 2 # Physical dimension for the 1D transverse and longitudinal field Ising model is 2

    # zero matrix, identity matrix and Pauli matrices X and Y - notice this implicitly assumes d = 2 for our degrees of freedom on a site
    zero = [0.0 0.0; 0.0 0.0]
    I = [1.0 0.0; 0.0 1.0]
    X = [0.0 1.0; 1.0 0.0]
    Z = [1.0 0.0; 0.0 -1.0]

    mpo = Vector{Array{ComplexF64}}(undef, N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor
    
    # mpo tensor at site 1
    
    mpo[1] = zeros((1, D, d, d))

    mpo[1][1,1,:,:] = g_x*X+g_z*Z
    mpo[1][1,2,:,:] = J*Z
    mpo[1][1,3,:,:] = I

    # mpo tensor at site N

    mpo[N] = zeros((D, 1, d, d))

    mpo[N][1,1,:,:] = I
    mpo[N][2,1,:,:] = Z
    mpo[N][3,1,:,:] = g_x*X+g_z*Z

    # mpo tensors at sites 2 to N-1

    for i in 2:N-1
        
        mpo[i] = zeros((D, D, d, d))
        
        mpo[i][1,1,:,:] = I
        mpo[i][1,2,:,:] = zero
        mpo[i][1,3,:,:] = zero
        mpo[i][2,1,:,:] = Z
        mpo[i][2,2,:,:] = zero
        mpo[i][2,3,:,:] = zero
        mpo[i][3,1,:,:] = g_x*X+g_z*Z
        mpo[i][3,2,:,:] = J*Z
        mpo[i][3,3,:,:] = I

    end

    return mpo # Vector{Array{ComplexF64, N} where N} and the array at each site is Array{ComplexF64, 4}

end

function get_identity_MPO(N::Int64, d::Int64)::Vector{Array{ComplexF64}}

    """
    Creates the identity MPO for a lattice of N sites with d degrees of freedom per site.

    Note 1: See notes for this function on my personal website on how to derived the MPO representation.

    Inputs:

    N = number of lattice sites (Integer)

    d = number of degrees of freedom per site - eg: d = 2 in the Ising model of spin-1/2 for spin up and spin down along z (Integer)

    Outputs:

    mpo = identity mpo as a vector of N elements and each element is a 4-tensor storing indices in the order shown schematically below

          sigma_i                                              3
             |                                                 |
    alpha -- Wi -- beta , which is stored in the order: 1 -- mpo[i] -- 2
             |                                                 |
         sigma_i_dash                                          4 


    In other words each tensor on each site has indices W_alpha,beta,sigma_i,sigma_i_dash
    """

    mpo = Vector{Array{ComplexF64}}(undef, N)

    on_each_site = reshape(Matrix{ComplexF64}(I, d, d), (1, 1, d, d))

    for i in 1:N
        mpo[i] = on_each_site
    end

    return mpo

end

function get_spin_half_MPO(N::Int64, measure_axis::String)::Vector{Array{ComplexF64}}

    """
    Returns the MPO as a Vector of 4-tensors that represents the operator of total spin 1/2 of the lattice

    O^j = sum_i (S^j_i) = sum_i (sigma^j_i/2) 
    
    where j = x,y,z depending on which axis we want to measure and i runs from 1 to N over lattice sites. 

    Note 1: See my personal website documentation to see how the MPO representation is derived.

    Inputs:

    N = number of lattice sites (Integer)

    measure_axis = which axis to measure the spin (String) takes values "x", "y" or "z"

    Output:

    mpo = vector where each element is a 4-tensor representing the operator as a W tensor
    """

    zero = [0.0 0.0; 0.0 0.0]
    identity = [1.0 0.0;0.0 1.0]

    if measure_axis == "x"
        operator = [0.0 0.5;0.5 0.0] # sigma_x/2 = pauli x operator divided by 2
    elseif measure_axis == "y"
        operator = [0.0 -0.5im;0.5im 0.0] # sigma_y/2 = pauli y operator divided by 2
    elseif measure_axis == "z"
        operator = [0.5 0.0;0.0 -0.5] # sigma_z/2 = pauli z operator divided by 2
    end

    mpo = Vector{Array{ComplexF64}}(undef, N)

    D = 2 # Bond dimension of MPO - this will be the dimension of the virtual/bond indices left and right of the tensor in a diagram
    d = 2 # Physical dimension of MPO - this will be the dimension of the physical indices above and below the tensor in a diagram

    # # Tensor at site 1

    mpo[1] = zeros(ComplexF64, 1,D,d,d)
    mpo[1][1,1,:,:] = operator
    mpo[1][1,2,:,:] = identity

    # # Tensor at site N

    mpo[N] = zeros(ComplexF64, D,1,d,d)
    mpo[N][1,1,:,:] = identity
    mpo[N][2,1,:,:] = operator

    # # Tensor at sites 2 to N-1

    tmp = zeros(ComplexF64, D,D,d,d)
    tmp[1,1,:,:] = identity
    tmp[1,2,:,:] = zero
    tmp[2,1,:,:] = operator
    tmp[2,2,:,:] = identity

    for i in 2:N-1
        mpo[i] = tmp
    end

    return mpo

end

function get_Schwinger_Wilson_MPO(N::Int64, l_0::Float64, x::Float64, lambda::Float64, m_g_ratio::Float64)::Vector{Array{ComplexF64}}

    """
    Builds the MPO for the 1+1 Schwinger model Hamiltonian using Wilson fermions, Gauss' law and a Jordan Wigner transformation to 
    put the Hamiltonian in full spin formulation. For the derivation see here: https://takisangelides.wixsite.com/personal/notes-sharing

    Inputs: 

    N = number of physical lattice sites (Integer) - in the spin formulation the lattice is spread to 2N sites

    l_0 = background electric field, theta/2pi (float) 

    x = 1/(a^2 * g^2) where a is lattice spacing and g is coupling constant (float)

    lambda = penalty term's lagrange multiplier (float)

    m_g_ratio = mass of fermion divided by coupling constant = m/g (float)

    Outputs:

    mpo = MPO of Schwinger model Hamiltonian using Wilson fermions (Vector of tensors of complex numbers)

    """

    A = -2*1im*(sqrt(x)*m_g_ratio + x)
    B = -2*1im*x
    C = (N-1)*l_0^2 + lambda*N/2 + N*(N-1)/4
    D = 7
    d = 2

    I = [1 0; 0 1]
    Z = [1 0; 0 -1]
    PLUS = [0 1; 0 0]
    MINUS = [0 0; 1 0]

    mpo = Vector{Array{ComplexF64}}(undef, 2*N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor

    for n in 1:2*N
    
        if n % 2 == 0
        
            if n == 2*N

                mpo[n] = zeros((D, 1, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][6,1,:,:] = 0.5*lambda.*Z
                mpo[n][7,1,:,:] = C/(2*N).*I

            else

                mpo[n] = zeros((D, D, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][6,1,:,:] = 0.5*(N - n/2 + lambda).*Z
                mpo[n][7,1,:,:] = C/(2*N).*I + l_0*(N-n/2).*Z
                mpo[n][5,2,:,:] = Z
                mpo[n][4,3,:,:] = Z
                mpo[n][6,6,:,:] = I
                mpo[n][7,6,:,:] = Z
                mpo[n][7,7,:,:] = I

            end
        else

            if n == 1
                
                mpo[n] = zeros((1, D, d, d))
                mpo[n][1,1,:,:] = C/(2*N).*I + l_0*(N - n/2 - 1/2).*Z
                mpo[n][1,2,:,:] = A.*PLUS
                mpo[n][1,3,:,:] = -A.*MINUS
                mpo[n][1,4,:,:] = B.*MINUS
                mpo[n][1,5,:,:] = -B.*PLUS
                mpo[n][1,6,:,:] = Z
                mpo[n][1,7,:,:] = I
            
            else

                mpo[n] = zeros((D, D, d, d))
                mpo[n][7,1,:,:] = C/(2*N).*I + l_0*(N - n/2 - 1/2).*Z
                mpo[n][7,2,:,:] = A.*PLUS
                mpo[n][7,3,:,:] = -A.*MINUS
                mpo[n][7,4,:,:] = B.*MINUS
                mpo[n][7,5,:,:] = -B.*PLUS
                mpo[n][7,6,:,:] = Z
                mpo[n][7,7,:,:] = I
                mpo[n][6,6,:,:] = I
                mpo[n][1,1,:,:] = I
                mpo[n][2,2,:,:] = Z
                mpo[n][3,3,:,:] = Z
                mpo[n][6,1,:,:] = 0.5*(N - n/2 - 1/2 + lambda).*Z

            end

        end
    
    end
    
    return mpo

end

function get_Schwinger_Wilson_MPO_Stefan(N::Int64, l_0::Float64, x::Float64, lambda::Float64, m_g_ratio::Float64)::Vector{Array{ComplexF64}}

    I = [1 0; 0 1]
    Z = [1 0; 0 -1]
    PLUS = [0 1; 0 0]
    MINUS = [0 0; 1 0]

    D = 7
    d = 2

    mpo = Vector{Array{ComplexF64}}(undef, 2*N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor

    for n in 1:2*N

        c_p = N - ceil(n/2)
    
        if n % 2 == 0
        
            if n == 2*N

                mpo[n] = zeros((D, 1, d, d))
                mpo[n][1,1,:,:] = ((1/(2*N))*((N-1)*l_0^2 + (1/4)*N*(N-1)+(1/2)*lambda*N)).*I
                mpo[n][2,1,:,:] = -2*1im*(m_g_ratio*sqrt(x) + x).*MINUS
                mpo[n][3,1,:,:] = 2*1im*x.*MINUS
                mpo[n][4,1,:,:] = 2*1im*(m_g_ratio*sqrt(x) + x).*PLUS
                mpo[n][5,1,:,:] = -2*1im*x.*PLUS
                mpo[n][6,1,:,:] = (1/2)*(lambda).*Z
                mpo[n][7,1,:,:] = I

            else

                mpo[n] = zeros((D, D, d, d))
                mpo[n][1,7,:,:] = ((1/(2*N))*((N-1)*l_0^2 + (1/4)*N*(N-1)+(1/2)*lambda*N)).*I + (l_0*c_p).*Z
                mpo[n][2,7,:,:] = -2*1im*(m_g_ratio*sqrt(x) + x).*MINUS
                mpo[n][3,7,:,:] = 2*1im*x.*MINUS
                mpo[n][4,7,:,:] = 2*1im*(m_g_ratio*sqrt(x) + x).*PLUS
                mpo[n][5,7,:,:] = -2*1im*x.*PLUS
                mpo[n][6,7,:,:] = (1/2)*(c_p + lambda).*Z
                mpo[n][7,7,:,:] = I
                mpo[n][1,1,:,:] = I
                mpo[n][2,2,:,:] = Z
                mpo[n][1,6,:,:] = Z
                mpo[n][4,4,:,:] = Z
                mpo[n][6,6,:,:] = I

            end
        else

            if n == 1
                
                mpo[n] = zeros((1, D, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][1,2,:,:] = PLUS
                mpo[n][1,4,:,:] = MINUS
                mpo[n][1,6,:,:] = Z
                mpo[n][1,7,:,:] = ((1/(2*N))*((N-1)*l_0^2 + (1/4)*N*(N-1)+(1/2)*lambda*N)).*I + (l_0*c_p).*Z
                
            
            else

                mpo[n] = zeros((D, D, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][1,2,:,:] = PLUS
                mpo[n][1,4,:,:] = MINUS
                mpo[n][1,6,:,:] = Z
                mpo[n][1,7,:,:] = ((1/(2*N))*((N-1)*l_0^2 + (1/4)*N*(N-1)+(1/2)*lambda*N)).*I + (l_0*c_p).*Z
                mpo[n][2,3,:,:] = Z
                mpo[n][4,5,:,:] = Z
                mpo[n][6,6,:,:] = I
                mpo[n][6,7,:,:] = (1/2)*(c_p + lambda).*Z
                mpo[n][7,7,:,:] = I

            end

        end
    
    end
    
    return mpo

end

function get_Schwinger_Wilson_general_r_MPO(N::Int64, l_0::Float64, x::Float64, lambda::Float64, m_g_ratio::Float64, r::Float64)::Vector{Array{ComplexF64}}

    """
    Builds the MPO for the 1+1 Schwinger model Hamiltonian using Wilson fermions with general r, Gauss' law and a Jordan Wigner transformation to 
    put the Hamiltonian in full spin formulation. 

    Inputs: 

    N = number of physical lattice sites (Integer) - in the spin formulation the lattice is spread to 2N sites

    l_0 = background electric field, theta/2pi (float) 

    x = 1/(a^2 * g^2) where a is lattice spacing and g is coupling constant (float)

    lambda = penalty term's lagrange multiplier (float)

    m_g_ratio = mass of fermion divided by coupling constant = m/g (float)

    Outputs:

    mpo = MPO of Schwinger model Hamiltonian using Wilson fermions (Vector of tensors of complex numbers)

    """

    A = (1/(2*N))*((N-1)*l_0^2 + lambda*N/2 + N*(N-1)/4)
    B = 1im*x*(r - 1)
    C = 2*1im*(sqrt(x)*m_g_ratio + x*r)
    D = 1im*x*(r + 1)
    
    D_bond = 7
    d = 2

    I = [1 0; 0 1]
    Z = [1 0; 0 -1]
    PLUS = [0 1; 0 0]
    MINUS = [0 0; 1 0]

    mpo = Vector{Array{ComplexF64}}(undef, 2*N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor

    for n in 1:2*N
    
        if n % 2 == 0
        
            if n == 2*N

                mpo[n] = zeros((D_bond, 1, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][7,1,:,:] = A.*I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][6,1,:,:] = 0.5*lambda.*Z

            else

                mpo[n] = zeros((D_bond, D_bond, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][7,7,:,:] = I
                mpo[n][7,1,:,:] = A.*I + l_0*(N-n/2).*Z
                mpo[n][6,6,:,:] = I
                mpo[n][7,3,:,:] = B.*MINUS
                mpo[n][2,1,:,:] = MINUS
                mpo[n][7,2,:,:] = -B.*PLUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][4,4,:,:] = Z
                mpo[n][5,5,:,:] = Z
                mpo[n][6,1,:,:] = 0.5*(N-n/2+lambda).*Z
                mpo[n][7,6,:,:] = Z 

            end
        else

            if n == 1
                
                mpo[n] = zeros((1, D_bond, d, d))
                mpo[n][1,7,:,:] = I
                mpo[n][1,1,:,:] = A.*I + l_0*(N-n/2-0.5).*Z
                mpo[n][1,3,:,:] = C.*MINUS
                mpo[n][1,5,:,:] = -D.*MINUS
                mpo[n][1,2,:,:] = -C.*PLUS
                mpo[n][1,4,:,:] = D.*PLUS 
                mpo[n][1,6,:,:] = Z
                
            
            else

                mpo[n] = zeros((D_bond, D_bond, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][7,7,:,:] = I
                mpo[n][7,1,:,:] = A.*I + l_0*(N-n/2-0.5).*Z
                mpo[n][6,6,:,:] = I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][7,3,:,:] = C.*MINUS
                mpo[n][7,5,:,:] = -D.*MINUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][7,2,:,:] = -C.*PLUS
                mpo[n][7,4,:,:] = D.*PLUS
                mpo[n][4,2,:,:] = Z
                mpo[n][5,3,:,:] = Z
                mpo[n][6,1,:,:] = 0.5*(N-n/2-0.5+lambda).*Z
                mpo[n][7,6,:,:] = Z
            
            end

        end
    
    end
    
    return mpo

end

function get_free_Wilson_general_r_MPO(N::Int64, x::Float64, m_g_ratio::Float64, r::Float64)::Vector{Array{ComplexF64}}

    A = (1/(2*N))*((N-1)*l_0^2 + lambda*N/2 + N*(N-1)/4)
    B = 1im*x*(r - 1)
    C = 2*1im*(sqrt(x)*m_g_ratio + x*r)
    D = 1im*x*(r + 1)
    
    D_bond = 6
    d = 2

    I = [1 0; 0 1]
    Z = [1 0; 0 -1]
    PLUS = [0 1; 0 0]
    MINUS = [0 0; 1 0]

    mpo = Vector{Array{ComplexF64}}(undef, 2*N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor

    for n in 1:2*N
    
        if n % 2 == 0
        
            if n == 2*N

                mpo[n] = zeros((D_bond, 1, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][3,1,:,:] = PLUS

            else

                mpo[n] = zeros((D_bond, D_bond, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][6,6,:,:] = I
                mpo[n][6,3,:,:] = B.*MINUS
                mpo[n][2,1,:,:] = MINUS
                mpo[n][6,2,:,:] = -B.*PLUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][4,4,:,:] = Z
                mpo[n][5,5,:,:] = Z 

            end
        else

            if n == 1
                
                mpo[n] = zeros((1, D_bond, d, d))
                mpo[n][1,6,:,:] = I
                mpo[n][1,3,:,:] = C.*MINUS
                mpo[n][1,5,:,:] = -D.*MINUS
                mpo[n][1,2,:,:] = -C.*PLUS
                mpo[n][1,4,:,:] = D.*PLUS
                
            
            else

                mpo[n] = zeros((D_bond, D_bond, d, d))
                mpo[n][1,1,:,:] = I
                mpo[n][6,6,:,:] = I
                mpo[n][2,1,:,:] = MINUS
                mpo[n][6,3,:,:] = C.*MINUS
                mpo[n][6,5,:,:] = -D.*MINUS
                mpo[n][3,1,:,:] = PLUS
                mpo[n][6,2,:,:] = -C.*PLUS
                mpo[n][6,4,:,:] = D.*PLUS
                mpo[n][4,2,:,:] = Z
                mpo[n][5,3,:,:] = Z
            
            end

        end
    
    end
    
    return mpo

end

function get_local_charge_MPO(N::Int64, site::Int64)::Vector{Array{ComplexF64}}

    """
    The charge operator for the Schwinger model using Wilson fermions at a particular spin site is returned as an MPO

    Inputs:

    N = number of physical lattice sites

    site = the site on which to act with the charge operator, on the extended spin lattice

    Output:

    mpo = the mpo for a single charge operator on a given site
    """

    mpo = Vector{Array{ComplexF64}}(undef, 2*N)
    d = 2
    D = 2
    Z = [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im]
    I = [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]

    for i in 1:2*N
    
        if i == 2*site-1 || i == 2*site
            if i == 1
                mpo[i] = zeros((1, D, d, d))
                mpo[i][1,1,:,:] = Z/2
                mpo[i][1,2,:,:] = I
            elseif i == 2*N
                mpo[i] = zeros((D, 1, d, d))
                mpo[i][1,1,:,:] = I
                mpo[i][2,1,:,:] = Z/2
            else
                mpo[i] = zeros((D, D, d, d))
                mpo[i][1,1,:,:] = I
                mpo[i][2,2,:,:] = I
                mpo[i][2,1,:,:] = Z/2
            end
        else
            if i == 1
                mpo[i] = zeros((1, D, d, d))
                mpo[i][1,2,:,:] = I
            elseif i == 2*N
                mpo[i] = zeros((D, 1, d, d))
                mpo[i][1,1,:,:] = I
            else
                mpo[i] = zeros((D, D, d, d))
                mpo[i][1,1,:,:] = I
                mpo[i][2,2,:,:] = I
            end
        end 
    
    end

    return mpo

end

function get_local_spin_MPO(N::Int64, site::Int64)::Vector{Array{ComplexF64}}

    """
    Gets the MPO for the Z operator acting at a particular spin site

    Inputs:

    N = number of physical lattice sites

    site = the site on which to act with the Z operator on the extended spin lattice

    Output:

    mpo = the mpo for a single Z operator on a given site
    """

    mpo = Vector{Array{ComplexF64}}(undef, 2*N)
    d = 2
    D = 2
    Z = [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im]
    I = [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]

    for i in 1:2*N
        if i == 2*N
            if i == site
                mpo[i] = zeros((D, 1, d, d))
                mpo[i][2,1,:,:] = Z
            else
                mpo[i] = zeros((D, 1, d, d))
                mpo[i][1,1,:,:] = I
            end
        elseif i == 1
            if i == site
                mpo[i] = zeros((1, D, d, d))
                mpo[i][1,1,:,:] = Z
            else
                mpo[i] = zeros((1, D, d, d))
                mpo[i][1,2,:,:] = I
            end
        else
            if i == site
                mpo[i] = zeros((D, D, d, d))
                mpo[i][2,1,:,:] = Z
            else
                mpo[i] = zeros((D, D, d, d))
                mpo[i][1,1,:,:] = I
                mpo[i][2,2,:,:] = I
            end
        end
    end

    return mpo


end

function get_penalty_term_MPO(N::Int64, lambda::Float64)::Vector{Array}

    """
    The penalty term in the dimensionless Hamiltonian W = 2/ag^2 H enforces the total charge of the ground state to be 0.
    We do this using a lagrange multiplier lambda and add to W the term lambda (sum from n = 1 to N Q_n)^2 = P_lambda. We
    have W + P_lambda = W'. So this function returns the mpo for P_lambda. Note for Wilson fermions is
    Q_n = psi dagger_n,2 psi_n,2 + psi dagger_n,1 psi_n,1 - 1.

    Inputs:
    
    N = number of lattice sites if the number of physical sites is M then N = 2M for our spin formulation

    lambda = lagrange multiplier (float)

    Output:

    mpo = the mpo for the penalty term which enforces zero charge on the g.s. (Vector of arrays)

    """

    mpo = Vector{Array}(undef, N) # An MPO is stored as a vector and each element of the vector stores a 4-tensor

    I = [1 0;0 1]
    Z = [1 0;0 -1]
    D = 3
    d = 2

    for n in 1:N
    
        if n == 1
            
            mpo[n] = zeros((1, D, d, d))
            mpo[n][1,1,:,:] = (lambda/4)*I
            mpo[n][1,2,:,:] = Z
            mpo[n][1,3,:,:] = I
            
        elseif n == N

            mpo[n] = zeros((D, 1, d, d))
            mpo[n][1,1,:,:] = I
            mpo[n][2,1,:,:] = (lambda/2)*Z
            mpo[n][3,1,:,:] = (lambda/4)*I 

        else

            mpo[n] = zeros((D, D, d, d))
            mpo[n][1,1,:,:] = I
            mpo[n][2,2,:,:] = I
            mpo[n][3,3,:,:] = I
            mpo[n][3,1,:,:] = (lambda/4)*I
            mpo[n][2,1,:,:] = (lambda/2)*Z
            mpo[n][3,2,:,:] = Z

        end

    end
    
    return mpo

end

function get_chiral_condensate_MPO(N::Int64)::Vector{Array{ComplexF64}}

    """
    Gets the MPO of the chiral condensate operator in the Wilson Schwinger model corresponding
    to the expression <psi_bar psi>.

    Input:

    N = number of spin sites meaning this is double the number of physical sites (Int)

    Output:

    mpo = the chiral condensate mpo (Vector of array of complex numbers)
    """

    mpo = Vector{Array{ComplexF64}}(undef, N)

    D = 4
    d = 2
    I = [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
    PLUS = [0.0+0.0im 1.0+0.0im; 0.0+0.0im 0.0+0.0im]
    MINUS = [0.0+0.0im 0.0+0.0im; 1.0+0.0im 0.0+0.0im]

    tensor_first = zeros(ComplexF64, 1, D, d, d)
    tensor_first[1,2,:,:] = -1im.*PLUS
    tensor_first[1,3,:,:] = 1im.*MINUS
    tensor_first[1,D,:,:] = I

    tensor_last = zeros(ComplexF64, D, 1, d, d)
    tensor_last[1,1,:,:] = I
    tensor_last[2,1,:,:] = MINUS
    tensor_last[3,1,:,:] = PLUS

    tensor_even = zeros(ComplexF64, D, D, d, d)
    tensor_even[1,1,:,:] = I
    tensor_even[2,1,:,:] = MINUS
    tensor_even[3,1,:,:] = PLUS
    tensor_even[D,D,:,:] = I

    tensor_odd = zeros(ComplexF64, D, D, d, d)
    tensor_odd[1,1,:,:] = I
    tensor_odd[D,2,:,:] = -1im.*PLUS
    tensor_odd[D,3,:,:] = 1im.*MINUS
    tensor_odd[D,D,:,:] = I

    for i in 1:N
        if i == 1
            mpo[i] = tensor_first
        elseif i == N
            mpo[i] = tensor_last
        elseif i % 2 == 0
            mpo[i] = tensor_even
        else
            mpo[i] = tensor_odd
        end
    end

    return mpo
end

function get_pseudo_momentum_MPO(N::Int64, x::Float64)::Vector{Array{ComplexF64}}

    """
    Inputs:

    N = number of physical lattice sites, as opposed to spin lattice sites which is 2N
    """

    mpo = Vector{Array{ComplexF64}}(undef, 2*N)

    D = 6
    d = 2
    I = [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
    PLUS = [0.0+0.0im 1.0+0.0im; 0.0+0.0im 0.0+0.0im]
    MINUS = [0.0+0.0im 0.0+0.0im; 1.0+0.0im 0.0+0.0im]
    Z = [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im]

    tensor_first = zeros(ComplexF64, 1, D, d, d)
    tensor_first[1,4,:,:] = MINUS
    tensor_first[1,5,:,:] = PLUS
    tensor_first[1,D,:,:] = I

    tensor_last = zeros(ComplexF64, D, 1, d, d)
    tensor_last[1,1,:,:] = I
    tensor_last[2,1,:,:] = 1im * x * MINUS
    tensor_last[3,1,:,:] = -1im * x * PLUS

    tensor = zeros(ComplexF64, D, D, d, d)
    tensor[1,1,:,:] = I
    tensor[2,1,:,:] = 1im * x * MINUS
    tensor[3,1,:,:] = -1im * x * PLUS
    tensor[3,1,:,:] = -1im * x * PLUS
    tensor[D,4,:,:] = MINUS
    tensor[D,5,:,:] = PLUS
    tensor[5,2,:,:] = Z
    tensor[4,3,:,:] = Z
    tensor[D,D,:,:] = I

    for i in 1:2*N
        if i == 1
            mpo[i] = tensor_first
        elseif i == 2*N
            mpo[i] = tensor_last
        else
            mpo[i] = tensor
        end
    end

    return mpo

end