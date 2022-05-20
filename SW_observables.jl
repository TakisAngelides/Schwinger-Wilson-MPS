include("MPO.jl")
include("variational_first_excited_state_MPS_algorithm.jl")

N = parse(Int, ARGS[1]) # Number of physical lattice sites
x = parse(Float64, ARGS[2]) # 1/(ag)^2
mg = parse(Float64, ARGS[3]) # m/g
D = parse(Int64, ARGS[4]) # Bond dimension of MPS to be computed
l_0 = parse(Float64, ARGS[5]) # l_0 = theta/(2*pi)
lambda = parse(Float64, ARGS[6]) # Lagrange multiplier to enforce total charge squared to be 0
acc = parse(Float64, ARGS[7]) # Tolerance for stopping condition of the variational algorithm
ms = parse(Int64, ARGS[8]) # Maximum number of sweeps of the variational algorithm
r = parse(Float64, ARGS[9])
choice = parse(Int64, ARGS[10])

mps = h5_to_mps(N, x, D, l_0, mg, ms, acc, lambda, r)
path = "/lustre/fs23/group/nic/tangelides/SW_Observables"

if choice == 1 # Energy

    mpo = get_Schwinger_Wilson_general_r_MPO(N, l_0, x, lambda, mg, r)
    energy = get_mpo_expectation_value(mps, mpo)
    
    text_file_name = "/Energy/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(energy)")
    end

elseif choice == 2 # Entropy

    entropy = entanglement_entropy(mps, N)
    
    text_file_name = "/Entropy/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(entropy)")
    end

elseif choice == 3 # Chiral Condensate

    cc_mpo = get_chiral_condensate_MPO(2*N)
    cc = real(get_mpo_expectation_value(mps, cc_mpo))
    
    text_file_name = "/Chiral Condensate/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(cc)")
    end
    
end
