using Oscar
Qx, s = polynomial_ring(QQ, :s);

# For the rest of the computations, we only need the algebraic number sqrt(-3), hence we can pass to a smaller number field
M, q3 = number_field(s^2 + 3, "q3");

#=B, z = polynomial_ring(M, :z => (1:8));

W_gen_1B = B(z[3]*z[7] + 2*z[4]*z[8]);
W_gen_2B = B(z[2]*z[4] + 2*q3*z[1]*z[7]);
W_gen_3B = B(z[2]*z[3] - 4*q3*z[1]*z[8]);
W_gen_4B = B(z[2]^3 + 12*q3*z[7]^2);
W_gen_5B = B(z[1]*z[2]^2 - 6*z[4]*z[7]);
W_gen_6B = B(z[1]^2*z[2] - q3*z[4]^2);

F_B = free_module(B, 1);
print(is_regular_sequence([W_gen_1B, W_gen_2B, W_gen_3B, W_gen_4B, W_gen_5B, W_gen_6B], F_B));=#

# In order to compute the blowup along the Weil divisor, we need to define a new polynomial ring in 9 variables, where the last variable t is used to model the homogeneous coordinates of the blowup
Bt, z, t = polynomial_ring(M, :z => (1:8), :t => (1:1));

# Embed the previously computed ideals (see invariants.jl) into the new polynomial ring Bt used to compute the blowups
rel_1 = q3*z[1]^3*z[5] - z[1]*z[3]*z[4] - 2*z[2]*z[6] - z[5]*z[8]; 
rel_2 = q3*z[1]^3*z[4] + z[1]*z[2]*z[5] - 2*z[3]*z[7] - z[4]*z[8];
rel_3 = z[1]*z[5]^2 + z[3]*z[8] + 2*z[4]*z[6];
rel_4 = z[1]*z[4]^2 - z[2]*z[8] - 2*z[5]*z[7];
rel_5 = -z[1]^4 - 3*q3*z[1]*z[8] + z[2]*z[3] - z[4]*z[5];
rel_6 = z[1]^3*z[8] + 1//4*z[2]*z[5]^2 - 1//4*z[3]*z[4]^2 - q3*z[6]*z[7] + q3*z[8]^2;
rel_7 = -2*z[1]^3*z[6] + q3*z[1]^2*z[3]*z[5] - z[3]^2*z[4] + z[5]^3 - 6*q3*z[6]*z[8];
rel_8 = -2*z[1]^3*z[7] + q3*z[1]^2*z[2]*z[4] + z[2]^2*z[5] - z[4]^3 - 6*q3*z[7]*z[8];
rel_9 = z[1]^3*z[8] + q3*z[1]^2*z[4]*z[5] + z[2]*z[5]^2 - z[3]*z[4]^2 + 3*q3*z[8]^2;
relations_2 = ideal(Bt, [rel_1, rel_2, rel_3, rel_4, rel_5, rel_6, rel_7, rel_8, rel_9]);

# Define the relations as given in the paper to check that they generate the same ideal
rel_pap_1 = q3*z[1]^3*z[5] - z[1]*z[3]*z[4] - 2*z[2]*z[6] - z[5]*z[8];
rel_pap_2 = q3*z[1]^3*z[4] + z[1]*z[2]*z[5] - 2*z[3]*z[7] - z[4]*z[8];
rel_pap_3 = z[1]*z[5]^2 + 2*z[4]*z[6] + z[3]*z[8];
rel_pap_4 = z[1]*z[4]^2 - 2*z[5]*z[7] - z[2]*z[8];
rel_pap_5 = -z[1]^4 + z[2]*z[3] - z[4]*z[5] - 3*q3*z[1]*z[8];
rel_pap_6 = z[1]^2*z[4]*z[5] + q3*z[1]^3*z[8] + 4*z[6]*z[7] - z[8]^2;
rel_pap_7 = q3*z[1]^2*z[3]*z[5] - 2*z[1]^3*z[6] - z[3]^2*z[4] + z[5]^3 - 6*q3*z[6]*z[8];
rel_pap_8 = q3*z[1]^2*z[2]*z[4] - 2*z[1]^3*z[7] - z[4]^3 + z[2]^2*z[5] - 6*q3*z[7]*z[8];
rel_pap_9 = 4*z[1]^2*z[4]*z[5] + q3*z[3]*z[4]^2 - q3*z[2]*z[5]^2 + 4*z[6]*z[7] + 8*z[8]^2;
relations = ideal(Bt, [rel_pap_1, rel_pap_2, rel_pap_3, rel_pap_4, rel_pap_5, rel_pap_6, rel_pap_7, rel_pap_8, rel_pap_9]);
print(relations == relations_2);  # returns true, hence both sets of relations generate the same ideal

W_gen_1 = z[3]*z[7] + 2*z[4]*z[8];
W_gen_2 = z[2]*z[4] + 2*q3*z[1]*z[7];
W_gen_3 = z[2]*z[3] - 4*q3*z[1]*z[8];
W_gen_4 = z[2]^3 + 12*q3*z[7]^2;
W_gen_5 = z[1]*z[2]^2 - 6*z[4]*z[7];
W_gen_6 = z[1]^2*z[2] - q3*z[4]^2;
weildivisor = ideal(Bt, [W_gen_1, W_gen_2, W_gen_3, W_gen_4, W_gen_5, W_gen_6]);

# Define the vector used to construct the underlying map of the blowup along the Weil divisor W_2
blowup_vector = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]] , gens(weildivisor*t[1]));

# Define the ring C=B[b_1,...,b_6] used to compute the rees algebra of the blowup along the Weil divisor W_2
C, z, b = polynomial_ring(M, :z => (1:8), :b => (1:6));
blowup_map = hom(C, Bt, blowup_vector);
kernel_blowup_map = preimage(blowup_map, relations);
basis_blowup = standard_basis(kernel_blowup_map, ordering=negdegrevlex(C));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal = ideal(C, collect(basis_blowup));

# Check that the blowup ideal is homogeneous with respect to the b-variables, i.e., in the ring B[b_1,...,b_6], so the z-variables have degree 0
# Note that this is done by evaluating all z-variables at 1 and checking homogeneity of the resulting polynomials in the b-variables because otherwise we would need to work 
# in a graded polynomial ring which does not allow the non-homogeneous generators needed to analyze the charts later on
b_part_kernel = [evaluate(gens(kernel_blowup_map)[i], [1, 2, 3, 4, 5, 6, 7, 8], [1, 1, 1, 1, 1, 1, 1, 1]) for i in eachindex(gens(kernel_blowup_map))];
check_homogeneity = [is_homogeneous(b_part_kernel[i]) for i in eachindex(b_part_kernel)];  # Check that all generators are homogeneous with respect to the b-variables
print(all(check_homogeneity));
b_part_kernel_basis = [evaluate(gens(basis_blowup_ideal)[i], [1, 2, 3, 4, 5, 6, 7, 8], [1, 1, 1, 1, 1, 1, 1, 1]) for i in eachindex(gens(basis_blowup_ideal))];
check_homogeneity_basis = [is_homogeneous(b_part_kernel_basis[i]) for i in eachindex(b_part_kernel_basis)];  # Check that all generators are homogeneous with respect to the b-variables
print(all(check_homogeneity_basis));

#=f_transform = hom(B, C, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]);
W_lift = [f_transform(W_gen_1),
          f_transform(W_gen_2),
          f_transform(W_gen_3),
          f_transform(W_gen_4),
          f_transform(W_gen_5),
          f_transform(W_gen_6)];

blowup_relations = [];
for i in 1:5
    for j in i+1:6
        push!(blowup_relations, W_lift[i]*b[j] - W_lift[j]*b[i])
    end
end
blowup_ideal = ideal(C, blowup_relations);=#

# Define the blowup ideal generated by the relations after computing a minimal generating set in Singular
kern_singular = ideal(C, [
    6*z[7]*b[2] - z[1]*b[4] + z[2]*b[5],
    3*z[4]*b[2] + q3*z[1]*b[5] - q3*z[2]*b[6],
    z[3]*b[2] + 2*z[4]*b[3] - q3*z[5]*b[6],
    z[5]*b[1] + z[6]*b[2] - z[8]*b[3],
    z[2]*b[1] - 2*z[8]*b[2] - z[7]*b[3],
    6*z[1]*b[1] + q3*z[3]*b[2] - q3*z[4]*b[3],
    z[1]*b[6]^2 + 2*b[1]*b[2] - b[3]*b[5],
    3*z[1]*z[4]*b[6] + 2*q3*z[2]*b[1] - q3*z[5]*b[5],
    z[1]^2*b[6] + z[5]*b[2] - z[2]*b[3],
    4*z[1]*z[4]*b[5] - z[2]*z[4]*b[6] + 2*q3*z[1]*z[7]*b[6] - q3*z[5]*b[4],
    4*z[5]*b[2]*b[3] - z[2]*b[3]^2 + z[3]*b[6]^2 + 4*q3*b[1]^2,
    z[1]*z[5]*b[3] + z[3]*b[1] + q3*z[6]*b[6],
    z[2]^2*b[3] + 12*q3*z[7]*b[1] - z[3]*b[4] + 4*q3*z[8]*b[5],
    z[1]*z[2]*b[3] - 6*z[4]*b[1] - z[3]*b[5] + 4*q3*z[8]*b[6],
    4*z[5]*b[2]^2 - z[2]*b[2]*b[3] + z[4]*b[6]^2 + 2*q3*b[1]*b[5],
    3*z[2]*b[2]^2 + q3*b[5]^2 - q3*b[4]*b[6],
    6*z[5]^2*b[2] - 3*z[2]*z[5]*b[3] - q3*z[3]*z[4]*b[6] + 4*q3*z[8]*b[1] - 2*q3*z[6]*b[5],
    3*z[2]*z[5]*b[2] - q3*z[4]^2*b[6] + 4*q3*z[7]*b[1] + 4*q3*z[8]*b[5],
    z[1]*z[5]*b[2] - 2*z[4]*b[1] + q3*z[8]*b[6],
    z[2]^2*b[2] - z[4]*b[4] - 2*q3*z[7]*b[5],
    z[1]*z[2]*b[2] - z[4]*b[5] - 2*q3*z[7]*b[6],
    z[1]*z[5]^2 + 2*z[4]*z[6] + z[3]*z[8],
    z[1]*z[4]^2 - 2*z[5]*z[7] - z[2]*z[8],
    2*z[1]^2*b[2]*b[5] - 2*z[1]*z[2]*b[2]*b[6] + q3*z[7]*b[6]^2 - q3*b[1]*b[4],
    z[1]^3*b[5] - z[1]^2*z[2]*b[6] - 4*q3*z[7]*b[1] - q3*z[8]*b[5],
    z[1]^3*b[4] - z[1]^2*z[2]*b[5] - 3*q3*z[4]^2*b[5] + 12*z[4]*z[7]*b[6] + 3*q3*z[8]*b[4],
    2*z[1]^2*z[7]*b[3] - 2*q3*z[4]^2*b[1] - q3*z[3]*z[4]*b[5] + z[3]*z[7]*b[6] - 6*z[4]*z[8]*b[6] + q3*z[6]*b[4],
    z[1]^2*z[4]*b[3] - q3*z[5]^2*b[2] - z[3]*z[4]*b[6] + 6*z[8]*b[1],
    z[1]^2*z[3]*b[3] - q3*z[5]^2*b[3] - z[3]^2*b[6] - 12*z[6]*b[1],
    z[1]^3*b[3] - z[1]*z[3]*b[6] + 2*q3*z[6]*b[2] + q3*z[8]*b[3],
    4*z[1]^2*b[2]^2 - 4*z[4]*b[2]*b[6] + q3*z[2]*b[6]^2 - q3*b[3]*b[4],
    z[1]^3*b[2] + 2*z[1]*z[4]*b[6] + 3*q3*z[8]*b[2] - q3*z[5]*b[5],
    4*z[1]^3*z[8] - z[3]*z[4]^2 + z[2]*z[5]^2 - 4*q3*z[6]*z[7] + 4*q3*z[8]^2,
    4*z[1]^2*z[4]*z[5] + q3*z[3]*z[4]^2 - q3*z[2]*z[5]^2 + 4*z[6]*z[7] + 8*z[8]^2,
    3*z[1]^2*z[3]*z[5] + 2*q3*z[1]^3*z[6] + q3*z[3]^2*z[4] - q3*z[5]^3 - 18*z[6]*z[8],
    3*z[1]^3*z[5] + q3*z[1]*z[3]*z[4] + 2*q3*z[2]*z[6] + q3*z[5]*z[8],
    3*z[1]^2*z[2]*z[4] + 2*q3*z[1]^3*z[7] + q3*z[4]^3 - q3*z[2]^2*z[5] - 18*z[7]*z[8],
    3*z[1]^3*z[4] - q3*z[1]*z[2]*z[5] + 2*q3*z[3]*z[7] + q3*z[4]*z[8],
    z[1]^4 - z[2]*z[3] + z[4]*z[5] + 3*q3*z[1]*z[8]
]);
print(kern_singular == basis_blowup_ideal);
b_part_kernel_sing = [evaluate(gens(kern_singular)[i], [1, 2, 3, 4, 5, 6, 7, 8], [1, 1, 1, 1, 1, 1, 1, 1]) for i in eachindex(gens(kern_singular))];
check_homogeneity_sing = [is_homogeneous(b_part_kernel_sing[i]) for i in eachindex(b_part_kernel)];  # Check that all generators are homogeneous with respect to the b-variables
print(all(check_homogeneity_sing));


# Function that computes a minimal embedding of a given ideal by eliminating variables that appear linearly
function minimal_embedding(I::MPolyIdeal)
    R = base_ring(I)
    # Work with a list of generators, that is modifidied during the process
    current_gens = collect(gens(I))
    
    # Store the remaining varibales in a list
    ring_vars = gens(R)
    active_vars = collect(ring_vars)
    
    # Save the eliminated variables and their substitution expressions in a dictionary
    eliminated_map = Dict{elem_type(R), elem_type(R)}()
    
    cleaning = true
    while cleaning
        cleaning = false  # Only set to true if we find a variable to eliminate
        
        # Iterate over the current generators to find a linear relation
        for g in current_gens
            if is_zero(g); continue; end
            
            found_v = nothing
            subst_expr = nothing
            
            # Check each active variable for linearity in the chosen generator g
            for v in active_vars
                if degree(g, v) == 1
                    c = derivative(g, v)  # Since g is linear in v, the derivative gives the coefficient (this is easier than extracting the coefficient as needed below)
                    
                    # We can only eliminate v if c is a non-zero constant
                    if is_constant(c) && !is_zero(c)
                        # Then we have c*v + rest = 0 => v = -rest/c and rest is everything in g except the term c*v
                        rest = g - c*v 
                        
                        subst_expr = -rest * inv(c)
                        found_v = v
                        break  # We found a variable to eliminate and thus can break the inner for loop
                    end
                end
            end
            
            if found_v !== nothing
                eliminated_map[found_v] = subst_expr
                filter!(x -> x != found_v, active_vars)
                
                # Search for the index of the variable that is set to be eliminated
                idx = findfirst(isequal(found_v), ring_vars)
                
                new_gens = elem_type(R)[]
                for poly in current_gens
                    push!(new_gens, evaluate(poly, [idx], [subst_expr]))
                end
                current_gens = new_gens
                
                print("Eliminating variable: $(found_v)\n")
                cleaning = true
                break  # Restart the while Loop
            end
        end
    end
    
    # Result: a triple of lists, consisting of the remaining variables, the ideal and the relations of the eliminations
    return active_vars, ideal(R, current_gens), eliminated_map
end



# Chart 1: b[1] = 1
kernel_1 = ideal(C, union(gens(kernel_blowup_map), [b[1] - 1]));
minimal_vars_1, minimal_ideal_1, elim_rules_1 = minimal_embedding(kernel_1);
print(minimal_vars_1);
kernel_basis_1 = ideal(C, union(gens(basis_blowup_ideal), [b[1] - 1]));
b_minimal_vars_1, b_minimal_ideal_1, b_elim_rules_1 = minimal_embedding(kernel_basis_1);
print(b_minimal_vars_1);


# Chart 2: b[2] = 1
kernel_2 = ideal(C, union(gens(kernel_blowup_map), [b[2] - 1]));
minimal_vars_2, minimal_ideal_2, elim_rules_2 = minimal_embedding(kernel_2);
print(minimal_vars_2);
kernel_basis_2 = ideal(C, union(gens(basis_blowup_ideal), [b[2] - 1]));
b_minimal_vars_2, b_minimal_ideal_2, b_elim_rules_2 = minimal_embedding(kernel_basis_2);
print(b_minimal_vars_2);


# Chart 3: b[3] = 1
kernel_3 = ideal(C, union(gens(kernel_blowup_map), [b[3] - 1]));
minimal_vars_3, minimal_ideal_3, elim_rules_3 = minimal_embedding(kernel_3);
print(minimal_vars_3);
kernel_basis_3 = ideal(C, union(gens(basis_blowup_ideal), [b[3] - 1]));
b_minimal_vars_3, b_minimal_ideal_3, b_elim_rules_3 = minimal_embedding(kernel_basis_3);
print(b_minimal_vars_3);


# Chart 4: b[4] = 1
kernel_4 = ideal(C, union(gens(kernel_blowup_map), [b[4] - 1]));
minimal_vars_4, minimal_ideal_4, elim_rules_4 = minimal_embedding(kernel_4);
print(minimal_vars_4);
kernel_basis_4 = ideal(C, union(gens(basis_blowup_ideal), [b[4] - 1]));
b_minimal_vars_4, b_minimal_ideal_4, b_elim_rules_4 = minimal_embedding(kernel_basis_4);
print(b_minimal_vars_4);


# Chart 5: b[5] = 1
kernel_5 = ideal(C, union(gens(kernel_blowup_map), [b[5] - 1]));
minimal_vars_5, minimal_ideal_5, elim_rules_5 = minimal_embedding(kernel_5);
print(minimal_vars_5);
kernel_basis_5 = ideal(C, union(gens(basis_blowup_ideal), [b[5] - 1]));
b_minimal_vars_5, b_minimal_ideal_5, b_elim_rules_5 = minimal_embedding(kernel_basis_5);
print(b_minimal_vars_5);


# Chart 6: b[6] = 1
kernel_6 = ideal(C, union(gens(kernel_blowup_map), [b[6] - 1]));
minimal_vars_6, minimal_ideal_6, elim_rules_6 = minimal_embedding(kernel_6);
print(minimal_vars_6);
kernel_basis_6 = ideal(C, union(gens(basis_blowup_ideal), [b[6] - 1]));
b_minimal_vars_6, b_minimal_ideal_6, b_elim_rules_6 = minimal_embedding(kernel_basis_6);
print(b_minimal_vars_6);


min_basis_1 = standard_basis(b_minimal_ideal_1, ordering=negdegrevlex(C))
min_basis_2 = standard_basis(b_minimal_ideal_2, ordering=negdegrevlex(C))
min_basis_3 = standard_basis(b_minimal_ideal_3, ordering=negdegrevlex(C))
min_basis_4 = standard_basis(b_minimal_ideal_4, ordering=negdegrevlex(C))
min_basis_5 = standard_basis(b_minimal_ideal_5, ordering=negdegrevlex(C))
min_basis_6 = standard_basis(b_minimal_ideal_6, ordering=negdegrevlex(C))


# Check the simplified ideals as given in the paper
c_2_test_ideal = ideal(C, (b[5]*b[6]-2*q3*z[1])^2+b[4]*(3*q3*b[3]-b[6]^3));
print(c_2_test_ideal == b_minimal_ideal_2);

c_3_test_ideal = ideal(C, [4*z[1]*b[1]+q3*z[3]*b[2]+z[5]*b[6], z[1]*z[5]+z[3]*b[1]+q3*z[6]*b[6], 
        z[1]^2*b[6]-z[3]*b[6]^2-4*q3*b[1]^2-3*z[5]*b[2], z[1]^2*z[3]-z[3]^2*b[6]-q3*z[5]^2-12*z[6]*b[1],
        z[1]^3-z[1]*z[3]*b[6]+q3*z[5]*b[1]+3*q3*z[6]*b[2]]);
print(c_3_test_ideal == b_minimal_ideal_3);

c_6_test_ideal = ideal(C, [z[4] + 2*q3*b[1]*b[5] + 4*z[5]*b[2]^2 + q3*z[4]*b[2]^2*b[3] + 2*b[1]*b[2]^2*b[3]*b[5] - b[2]*b[3]^2*b[5]^2,
        3*z[5] + 2*q3*z[4]*b[3] + 12*b[1]^2*b[2] - 4*q3*z[5]*b[2]^2*b[3] + 3*z[4]*b[2]^2*b[3]^2 - 2*q3*b[1]*b[2]^2*b[3]^2*b[5] + q3*b[2]*b[3]^3*b[5]^2]);
print(c_6_test_ideal == b_minimal_ideal_6);


function evaluate_variables(gen_list, idx)  # Note that the indices of b[1],...,b[6] are 9,...,14 in C
    zeros_vec = zeros(Int, length(idx))
    return [evaluate(gen_list[i], idx, zeros_vec) for i in eachindex(gen_list)]
end

function nonzero_constant_indices(gen_list)  # Returns the indices of non-zero constant generators, used to detect constant generators after evaluation
    return [i for i in eachindex(gen_list)
            if is_constant(gen_list[i]) && gen_list[i] != 0]
end


# Analyze singularities of the blowup charts
# Chart 1
print(b_minimal_ideal_1 == minimal_ideal_1);  # Returns false, so the embedding depends on the chosen basis

gen_1 = gens(minimal_ideal_1)[[3; 8; 13; 18:19]];
test_ideal_1 = ideal(C, gen_1);
print(minimal_ideal_1 == test_ideal_1);  # Shows that the ideal can be simplified to only 5 generators
eval_1 = evaluate_variables(gen_1, [11, 14]);
nonzero_eval_gen_1 = nonzero_constant_indices(eval_1);  # Shows that there is a non-zero constant generator when both b[3] and b[6] are set to zero
print(nonzero_eval_gen_1);

eval_b_ideal_1 = evaluate_variables(gens(b_minimal_ideal_1), [11, 14]);
nonzero_eval_b_ideal_1 = nonzero_constant_indices(eval_b_ideal_1);  # Shows that there is a non-zero constant generator when both b[3] and b[6] are set to zero
print(nonzero_eval_b_ideal_1);

b_gen_1 = gens(b_minimal_ideal_1)[[2; 7; 12; 13; 16; 18]];
test_b_ideal_1 = ideal(C, b_gen_1);
print(b_minimal_ideal_1 == test_b_ideal_1);  # Shows that the ideal can be simplified to only 6 generators
eval_b_min_1 = evaluate_variables(b_gen_1, [11, 14]);
nonzero_eval_b_min_1 = nonzero_constant_indices(eval_b_min_1);  # Shows that there is a non-zero constant generator when both b[3] and b[6] are set to zero
print(nonzero_eval_b_min_1);

#=
F_1 = free_resolution(minimal_ideal_1, algorithm = :mres);  # Runs out of time, but is not needed for the analysis

F_1_b = free_resolution(b_minimal_ideal_1, algorithm = :mres);  # Runs out of time, but is not needed for the analysis

X_1 = spec(C, minimal_ideal_1);
sing_1 = singular_locus(X_1);

X_1_b = spec(C, b_minimal_ideal_1);
sing_1_b = singular_locus(X_1_b);
=#


# Chart 2: calculate singular locus based on simplified ideal
print(b_minimal_ideal_2 == minimal_ideal_2);  # Returns true, so the embedding does not depend on the chosen basis
#=
F_2 = free_resolution(b_minimal_ideal_2, algorithm = :mres);
complex_2 = augmented_complex(F_2);
map_2 = map(complex_2, 0);
mat_2 = matrix(map_2);  # only one non-zero entry at index 33
generator_2 = gens(b_minimal_ideal_2)[33];
=#
generator_2 = 12*z[1]^2 + 4*q3*z[1]*b[5]*b[6] - 3*q3*b[3]*b[4] + b[4]*b[6]^3 - b[5]^2*b[6]^2;  # Simplified generator found by inspecting the standard basis
c_2_simplified_ideal = ideal(C, generator_2);
print(b_minimal_ideal_2 == c_2_simplified_ideal);  # Shows that the ideal can be simplified to only one generator

X_2 = spec(C, c_2_simplified_ideal);
sing_2 = singular_locus(X_2);
y_1 =  -3*q3*b[3] + b[6]^3;
y_2 =  -2*q3*z[1] + b[5]*b[6];
test_gen = b[4]*y_1 - y_2^2;  # Coordinate change to see the A1 singularity
print(test_gen == generator_2);  # Shows that the second chart has an A1 singularity after this coordinate change


# Chart 3: calculate singular locus based on simplified ideal
print(b_minimal_ideal_3 == minimal_ideal_3);  # Returns true, so the embedding does not depend on the chosen basis
F_3 = free_resolution(b_minimal_ideal_3, algorithm = :mres);
complex_3 = augmented_complex(F_3);
map_3 = map(complex_3, 0);
mat_3 = matrix(map_3);  # only non-zero entries at indices (1,12), (2,15), (3,19), (4,21), (5,22)
generator_3_1 = gens(b_minimal_ideal_3)[12];
generator_3_2 = gens(b_minimal_ideal_3)[15];
generator_3_3 = gens(b_minimal_ideal_3)[19];
generator_3_4 = gens(b_minimal_ideal_3)[21];
generator_3_5 = gens(b_minimal_ideal_3)[22];
c_3_simplified_ideal = ideal(C, [generator_3_1, generator_3_2, generator_3_3, generator_3_4, generator_3_5]);
print(b_minimal_ideal_3 == c_3_simplified_ideal);  # Shows that the ideal can be simplified to only five generators
#X_3 = spec(C, c_3_simplified_ideal);
#sing_3 = singular_locus(X_3);
# We want to show that this chart is isomorphic to the four dimensional S_3-quotient using a description established by Weyl (using polarizations of elementary symmetric polynomials)
R_S3, x, y = polynomial_ring(M, :x => 1:3, :y => 1:3);
p_1 = x[1]*x[2] + x[1]*x[3] + x[2]*x[3];
p_2 = x[1]*y[2] + x[1]*y[3] + x[2]*y[1] + x[2]*y[3] + x[3]*y[1] + x[3]*y[2];
p_3 = y[1]*y[2] + y[1]*y[3] + y[2]*y[3];
p_4 = x[1]*x[2]*x[3];
p_5 = x[1]*x[2]*y[3] + x[1]*x[3]*y[2] + x[2]*x[3]*y[1];
p_6 = x[1]*y[2]*y[3] + x[2]*y[1]*y[3] + x[3]*y[1]*y[2];
p_7 = y[1]*y[2]*y[3];
ker_1 = x[1] + x[2] + x[3];
ker_2 = y[1] + y[2] + y[3];
S3_kernel = ideal(R_S3, [ker_1, ker_2]);
S3_invariants_list = [p_1, p_2, p_3, p_4, p_5, p_6, p_7];
# Using these invariants, we compute the relations among them to get a description of the S_3-quotient in C^7
S_S3, u = polynomial_ring(M, :u => 1:7);
f_S3 = hom(S_S3, R_S3, S3_invariants_list);
S3_relations = preimage(f_S3, S3_kernel);
standard_basis_S3 = standard_basis(S3_relations, ordering=negdegrevlex(S_S3));
# Inspection of the standard basis shows that the relations can be generated by the following five polynomials:
t_1 = 3*u[1]*u[7] - u[2]*u[6] + u[3]*u[5];
t_2 = u[1]*u[6] - u[2]*u[5] + 3*u[3]*u[4];
t_3 = -4*u[1]*u[3]^2 + u[2]^2*u[3] + 9*u[5]*u[7] - 3*u[6]^2;
t_4 = -4*u[1]*u[2]*u[3] + u[2]^3 + 27*u[4]*u[7] - 3*u[5]*u[6];
t_5 = -4*u[1]^2*u[3] + u[1]*u[2]^2 + 9*u[4]*u[6] - 3*u[5]^2;
S3_simplified_ideal = ideal(S_S3, [t_1, t_2, t_3, t_4, t_5]);
print(S3_relations == S3_simplified_ideal);
# Now define a map describing an isomorphism between the third chart of the blowup and the S_3-quotient by a linear change of coordinates
# Note that the eliminated variables are sent to zero
f_blowup_S3 = hom(C, S_S3, [u[2], 0, (-1//q3)*u[3], 0, -u[6], 1//4*u[7], 0, 0, -q3*u[5], -12*q3*u[4], 0, 0, 0, -4*q3*u[1]]);
image_chart3 = f_blowup_S3(c_3_simplified_ideal);
print(image_chart3 == S3_simplified_ideal);  # Shows that the third chart is isomorphic to the four dimensional S_3-quotient


# Chart 4
print(b_minimal_ideal_4 == minimal_ideal_4);  # Returns true, so the embedding does not depend on the chosen basis
X_4 = spec(C, b_minimal_ideal_4);
sing_4 = singular_locus(X_4);
print(is_smooth(X_4));  # Returns true


# Chart 5
print(b_minimal_ideal_5 == minimal_ideal_5);  # Returns false, so the embedding depends on the chosen basis
rad_5 = radical(minimal_ideal_5);
eval_rad_5 = evaluate_variables(gens(rad_5), [10, 12]);
nonzero_eval_rad_5 = nonzero_constant_indices(eval_rad_5);  # Shows that there is a non-zero constant generator when both b[2] and b[4] are
print(nonzero_eval_rad_5);

rad_b_5 = radical(b_minimal_ideal_5);
eval_rad_b_5 = evaluate_variables(gens(rad_b_5), [10, 12]);
nonzero_eval_rad_b_5 = nonzero_constant_indices(eval_rad_b_5);  # Shows that there is a non-zero constant generator when both b[2] and b[4] are set to zero
print(nonzero_eval_rad_b_5);
#=
F_5 = free_resolution(b_minimal_ideal_5, algorithm = :mres);
complex_5 = augmented_complex(F_5);
map_5 = map(complex_5, 0);
X_5 = spec(C, b_minimal_ideal_5);
sing_5 = singular_locus(X_5);
=#


# Chart 6: calculate singular locus based on simplified ideal
print(b_minimal_ideal_6 == minimal_ideal_6);  # Returns false, so the embedding depends on the chosen basis
#=
F_6 = free_resolution(b_minimal_ideal_6, algorithm = :mres);
complex_6 = augmented_complex(F_6);
map_6 = map(complex_6, 0);
mat_6 = matrix(map_6);  # only non-zero entries at indices (1,19), (2,38)
generator_6_1 = gens(b_minimal_ideal_6)[19];
generator_6_2 = gens(b_minimal_ideal_6)[38];
=#
gen_6_1 = z[4] + 2*q3*b[1]*b[5] + 4*z[5]*b[2]^2 + q3*z[4]*b[2]^2*b[3] + 2*b[1]*b[2]^2*b[3]*b[5] - b[2]*b[3]^2*b[5]^2;
gen_6_2 =  3*z[5] + 2*q3*z[4]*b[3] + 12*b[1]^2*b[2] - 4*q3*z[5]*b[2]^2*b[3] + 3*z[4]*b[2]^2*b[3]^2 - 2*q3*b[1]*b[2]^2*b[3]^2*b[5] + q3*b[2]*b[3]^3*b[5]^2;
c_6_simplified_ideal = ideal(C, [gen_6_1, gen_6_2]);  # Simplified generators found by inspecting the standard basis
print(b_minimal_ideal_6 == c_6_simplified_ideal);  # Shows that the ideal can be simplified to only two generators

X_6 = spec(C, c_6_simplified_ideal);
sing_6 = singular_locus(X_6);
sing_6_ideal = modulus(OO(sing_6[1]));
eval_sing_6 = evaluate_variables(gens(sing_6_ideal), [10, 11]);
nonzero_eval_sing_6 = nonzero_constant_indices(eval_sing_6);  # Shows that there is a non-zero constant generator when both b[2] and b[3] are set to zero
print(nonzero_eval_sing_6);


# Finally, we compute the fiber of the blowup map over the origin
#fiber_origin_ideal = ideal(C, union(gens(kernel_blowup_map), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
#reduced_fiber_origin_ideal = radical(fiber_origin_ideal);

fiber_origin_basis_ideal = ideal(C, union(gens(basis_blowup_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
reduced_fiber_origin_basis_ideal = radical(fiber_origin_basis_ideal);

X_fiber = spec(C, reduced_fiber_origin_basis_ideal);
components_fiber = irreducible_components(X_fiber);

X_nonred_fiber = spec(C, fiber_origin_basis_ideal);

min_basis_fiber = standard_basis(modulus(OO(X_nonred_fiber)), ordering=negdegrevlex(C));
X_min_basis_fiber = spec(C, ideal(C, min_basis_fiber));

sing_fiber = singular_locus(X_min_basis_fiber);
sing_ideal = modulus(OO(sing_fiber[1]));
reduced_sing_ideal = radical(sing_ideal);

min_basis_fiber_W_1 = [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], 12*b[1]^2 + 4*q3*z[4]*b[2]*b[3] - q3*z[3]*b[3]^2 + q3*z[2]*b[6]^2, 2*b[1]*b[2] + b[3]*b[5] - z[1]*b[6]^2,
    3*b[1]*b[4] - 3*z[6]*b[6]^2 + 2*q3*z[1]^2*b[2]*b[5] - 2*q3*z[1]*z[3]*b[2]*b[6], 3*b[3]*b[4] + 4*q3*z[5]*b[2]*b[6] - 3*z[3]*b[6]^2 - 4*q3*z[1]^2*b[2]^2,
    6*b[1]*b[5] - 4*q3*z[4]*b[2]^2 + q3*z[3]*b[2]*b[3] - q3*z[5]*b[6]^2, b[5]^2 - b[4]*b[6] + q3*z[3]*b[2]^2];
X_min_basis_fiber_W_1 = spec(C, ideal(C, min_basis_fiber_W_1));

min_basis_fiber_W_2 = [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], 12*b[1]^2 - 4*q3*z[5]*b[2]*b[3] + q3*z[2]*b[3]^2 - q3*z[3]*b[6]^2, 2*b[1]*b[2] - b[3]*b[5] + z[1]*b[6]^2, 
    3*b[1]*b[4] - 3*z[7]*b[6]^2 + 2*q3*z[1]^2*b[2]*b[5] - 2*q3*z[1]*z[2]*b[2]*b[6], 3*b[3]*b[4] - 4*q3*z[4]*b[2]*b[6] - 3*z[2]*b[6]^2 + 4*q3*z[1]^2*b[2]^2, 
    6*b[1]*b[5] - 4*q3*z[5]*b[2]^2 + q3*z[2]*b[2]*b[3] - q3*z[4]*b[6]^2, b[5]^2 - b[4]*b[6] - q3*z[2]*b[2]^2];
X_min_basis_fiber_W_2 = spec(C, ideal(C, min_basis_fiber_W_2));

common_fiber = intersect(X_min_basis_fiber_W_1, X_min_basis_fiber);
print(common_fiber == X_min_basis_fiber);
common_fiber_ideal = modulus(OO(common_fiber));
common_fiber_min_basis = standard_basis(common_fiber_ideal, ordering=negdegrevlex(C));

embed_B = hom(Bt, C, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], 0]);
embedded_relations = ideal(C, [embed_B(g) for g in gens(relations)]);
C_quo, pi_C = quo(C, embedded_relations);
quo_fiber_W_1 = ideal(C_quo, [pi_C(g) for g in min_basis_fiber_W_1]);
quo_fiber_W_2 = ideal(C_quo, [pi_C(g) for g in min_basis_fiber_W_2]);
print(quo_fiber_W_1 == quo_fiber_W_2);