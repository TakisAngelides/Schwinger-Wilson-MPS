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
D_p = parse(Int64, ARGS[9]) # Bond dimension of ansatz
mg_p = parse(Float64, ARGS[10]) # m/g of ansatz
r = parse(Float64, ARGS[11]) # Wilson parameter

mpo = get_Schwinger_Wilson_general_r_MPO(N, l_0, x, lambda, mg, r)
variational_ground_state_algorithm(2*N, 2, x, l_0, lambda, mg, ms, acc, D, D_p, mg_p, mpo, r) # Computes and saves an MPS