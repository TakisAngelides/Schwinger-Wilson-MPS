These set of files are a Julia implementation of mainly the variational ground state search algorithm using matrix product states. There is also an implementation to get
the first excited state. The reason the title of this repository is Schwinger Wilson MPS is because it implements the Schwinger model with Wilson fermions MPO and performs
the variational ground state search algorithm on it.

**Functions in variational_MPS_algorithm.jl:**

initialize_MPS

contraction

gauge_site

gauge_mps!

initialize_L_R_states

get_Heff

get_updated_site

update_states!

variational_ground_state_MPS

variational_ground_state_MPS_from_previous

**Functions in MPO.jl:**

get_Ising_MPO

get_identity_MPO

get_spin_half_MPO

get_Schwinger_Wilson_MPO

get_local_charge_MPO

get_penalty_term_MPO

**Functions in variational_first_excited_state_MPS_algorithm.jl:**

initialize_left_right_states

update_p_states!

get_H_projection

get_updated_site_first_excited

variational_first_excited_MPS

**Functions in utility_functions.jl:**

get_GHZ_mps

get_electric_field_configuration

operator_to_sparse_matrix

get_Schwinger_hamiltonian_matrix

get_penalty_term_matrix

mpo_to_matrix

act_mpo_on_mps

get_spin_half_expectation_value

get_mpo_expectation_value

inner_product_MPS

quantum_state_coefficients

gauge_site_center_orthogonal

center_orthogonalize!

norm_center_orthogonal

entanglement_entropy_inefficient

entanglement_entropy_old

entanglement_entropy

generate_Schwinger_data

generate_entropy_data
