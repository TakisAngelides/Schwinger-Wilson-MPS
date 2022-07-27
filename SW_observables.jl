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

elseif choice == 4 # Electric field 

    ef = real(sum(get_electric_field_configuration(l_0, mps)))
    
    text_file_name = "/Electric Field/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(ef)")
    end

elseif choice == 5 # Electric field but avoiding the lattice edges

    left_edge = floor(Int, N*0.48)
    right_edge = floor(Int, N*0.52)
    efl = get_electric_field_configuration(l_0, mps) # Electric field list
    middle_efl = efl[left_edge:right_edge] # Both are inclusive in Julia and if its [5:5] it will still take the element 5
    number_of_links = length(middle_efl)
    ef = real(sum(middle_efl))
    
    text_file_name = "/Electric Field/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(ef),$(number_of_links)")
    end

elseif choice == 6 # Particle number

    mpo_particle_number = get_particle_number_MPO(N)
    pn = real(get_mpo_expectation_value(mps, mpo_particle_number))
    
    text_file_name = "/Particle Number/N_$(N)_x_$(x)_D_$(D)_l0_$(l_0)_mg_$(mg)_ms_$(ms)_acc_$(acc)_lam_$(lambda)_r_$(r).txt"
    path_to_text_file = path*text_file_name
    
    open(path_to_text_file, "w") do f
        write(f, "$(pn)")
    end

end
