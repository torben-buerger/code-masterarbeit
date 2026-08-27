using Oscar
Qx, s = polynomial_ring(QQ, :s);

# For the rest of the computations, we only need the algebraic number sqrt(-3), hence we can pass to a smaller number field
M, q3 = number_field(s^2 + 3, "q3");

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

W_gen_1 = z[2]*z[6] + 2*z[5]*z[8];
W_gen_2 = z[3]*z[5] + 2*q3*z[1]*z[6];
W_gen_3 = z[2]*z[3] - 4*q3*z[1]*z[8];
W_gen_4 = z[3]^3 - 12*q3*z[6]^2;
W_gen_5 = z[1]*z[3]^2 + 6*z[5]*z[6];
W_gen_6 = z[1]^2*z[3] + q3*z[5]^2;
weildivisor = ideal(Bt, [W_gen_1, W_gen_2, W_gen_3, W_gen_4, W_gen_5, W_gen_6]);

# Define the vector used to construct the underlying map of the blowup along the Weil divisor W_1
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

eval_1 = evaluate_variables(gens(minimal_ideal_1), [11, 14]);
nonzero_eval_gen_1 = nonzero_constant_indices(eval_1);  # Shows that there is a non-zero constant generator when both b[3] and b[6] are set to zero
print(nonzero_eval_gen_1);

eval_b_ideal_1 = evaluate_variables(gens(b_minimal_ideal_1), [11, 14]);
nonzero_eval_b_ideal_1 = nonzero_constant_indices(eval_b_ideal_1);  # Shows that there is a non-zero constant generator when both b[3] and b[6] are set to zero
print(nonzero_eval_b_ideal_1);


# Chart 2: calculate singular locus based on simplified ideal
print(b_minimal_ideal_2 == minimal_ideal_2);  # Returns true, so the embedding does not depend on the chosen basis

generator_2 = 12*z[1]^2 + 3*q3*b[3]*b[4] - 4*q3*z[1]*b[5]*b[6] - b[5]^2*b[6]^2 + b[4]*b[6]^3;  # Simplified generator found by inspecting the standard basis
c_2_simplified_ideal = ideal(C, generator_2);
print(b_minimal_ideal_2 == c_2_simplified_ideal);  # Shows that the ideal can be simplified to only one generator

X_2 = spec(C, c_2_simplified_ideal);
sing_2 = singular_locus(X_2);
y_1 =  3*q3*b[3] + b[6]^3;
y_2 =  2*q3*z[1] + b[5]*b[6];
test_gen = b[4]*y_1 - y_2^2;  # Coordinate change to see the A1 singularity
print(test_gen == generator_2);  # Shows that the second chart has an A1 singularity after this coordinate change


# Chart 3: calculate singular locus based on simplified ideal
print(b_minimal_ideal_3 == minimal_ideal_3);  # Returns true, so the embedding does not depend on the chosen basis
F_3 = free_resolution(b_minimal_ideal_3, algorithm = :mres);
complex_3 = augmented_complex(F_3);
map_3 = map(complex_3, 0);
mat_3 = matrix(map_3);  # only non-zero entries at indices (1,11), (2,16), (3,19), (4,20), (5,23)
generator_3_1 = gens(b_minimal_ideal_3)[11];
generator_3_2 = gens(b_minimal_ideal_3)[16];
generator_3_3 = gens(b_minimal_ideal_3)[19];
generator_3_4 = gens(b_minimal_ideal_3)[20];
generator_3_5 = gens(b_minimal_ideal_3)[23];
c_3_simplified_ideal = ideal(C, [generator_3_1, generator_3_2, generator_3_3, generator_3_4, generator_3_5]);
print(b_minimal_ideal_3 == c_3_simplified_ideal);  # Shows that the ideal can be simplified to only five generators
X_3 = spec(C, c_3_simplified_ideal);
sing_3 = singular_locus(X_3);

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
f_blowup_S3 = hom(C, S_S3, [u[2], (-1//3)*q3*u[3], 0, u[6], 0, 0, 1//4*u[7], 0, q3*u[5], -12*q3*u[4], 0, 0, 0, 4*q3*u[1]]);
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


# Chart 6: calculate singular locus based on simplified ideal
print(b_minimal_ideal_6 == minimal_ideal_6);  # Returns false, so the embedding depends on the chosen basis

gen_6_1 = z[5] + 2*q3*b[1]*b[5] + 4*z[4]*b[2]^2 - q3*z[5]*b[2]^2*b[3] - 2*b[1]*b[2]^2*b[3]*b[5] - b[2]*b[3]^2*b[5]^2;
gen_6_2 =  3*z[4] - 2*q3*z[5]*b[3] + 12*b[1]^2*b[2] + 4*q3*z[4]*b[2]^2*b[3] + 3*z[5]*b[2]^2*b[3]^2 - 2*q3*b[1]*b[2]^2*b[3]^2*b[5] - q3*b[2]*b[3]^3*b[5]^2;
c_6_simplified_ideal = ideal(C, [gen_6_1, gen_6_2]);  # Simplified generators found by inspecting the standard basis
print(b_minimal_ideal_6 == c_6_simplified_ideal);  # Shows that the ideal can be simplified to only two generators

X_6 = spec(C, c_6_simplified_ideal);
sing_6 = singular_locus(X_6);
sing_6_ideal = modulus(OO(sing_6[1]));
eval_sing_6 = evaluate_variables(gens(sing_6_ideal), [10, 11]);
nonzero_eval_sing_6 = nonzero_constant_indices(eval_sing_6);  # Shows that there is a non-zero constant generator when both b[2] and b[3] are set to zero
print(nonzero_eval_sing_6);


# The charts admitting singularities are blown up in their respective reduced singular locus to obtain the resolution as described in the paper

sing_2_ideal = modulus(OO(sing_2[1]));
sing_3_ideal = modulus(OO(sing_3[1]));
sing_6_ideal = modulus(OO(sing_6[1]));

min_sing_2 = standard_basis(sing_2_ideal, ordering=negdegrevlex(C));
min_sing_3 = standard_basis(sing_3_ideal, ordering=negdegrevlex(C));
min_sing_6 = gens(radical(sing_6_ideal));

Ct, z, b, t = polynomial_ring(M, :z => (1:8), :b => (1:6), :t => (1:1));
embed_Ct = hom(C, Ct, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]]);

blowup_vec_2 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_Ct(f)*t[1] for f in min_sing_2]);
blowup_vec_3 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_Ct(f)*t[1] for f in min_sing_3]);
blowup_vec_6 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_Ct(f)*t[1] for f in min_sing_6]);

# We compute the blowups locally above the covering charts having a singularity by using the Rees algebra as before
D2, z, b, w = polynomial_ring(M, :z => (1:8), :b => (1:6), :w => (1:3));
blowup_map_2 = hom(D2, Ct, blowup_vec_2);
relations_chart_2 = ideal(Ct, [embed_Ct(f) for f in gens(c_2_simplified_ideal)]);
kernel_blowup_map_2 = preimage(blowup_map_2, relations_chart_2);
basis_blowup_2 = standard_basis(kernel_blowup_map_2, ordering=negdegrevlex(D2));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_2 = ideal(D2, collect(basis_blowup_2));
print(basis_blowup_ideal_2 == kernel_blowup_map_2);

# Chart 1: w[1] = 1
kernel_basis_21 = ideal(D2, union(basis_blowup_2, [w[1] - 1]));
b_minimal_vars_21, b_minimal_ideal_21, b_elim_rules_21 = minimal_embedding(kernel_basis_21);
X_21 = spec(D2, b_minimal_ideal_21);
print(is_smooth(X_21));  # returns true

# Chart 2: w[2] = 1
kernel_basis_22 = ideal(D2, union(basis_blowup_2, [w[2] - 1]));
b_minimal_vars_22, b_minimal_ideal_22, b_elim_rules_22 = minimal_embedding(kernel_basis_22);  # becomes the 0 ideal

# Chart 3: w[3] = 1
kernel_basis_23 = ideal(D2, union(basis_blowup_2, [w[3] - 1]));
b_minimal_vars_23, b_minimal_ideal_23, b_elim_rules_23 = minimal_embedding(kernel_basis_23);  # becomes the 0 ideal


D3, z, b, u = polynomial_ring(M, :z => (1:8), :b => (1:6), :u => (1:17));
blowup_map_3 = hom(D3, Ct, blowup_vec_3);
relations_chart_3 = ideal(Ct, [embed_Ct(f) for f in gens(c_3_simplified_ideal)]);
kernel_blowup_map_3 = preimage(blowup_map_3, relations_chart_3);
basis_blowup_3 = standard_basis(kernel_blowup_map_3, ordering=negdegrevlex(D3));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_3 = ideal(D3, collect(basis_blowup_3));
print(basis_blowup_ideal_3 == kernel_blowup_map_3);

# Chart 1: u[1] = 1
kernel_basis_31 = ideal(D3, union(basis_blowup_3, [u[1] - 1]));
b_minimal_vars_31, b_minimal_ideal_31, b_elim_rules_31 = minimal_embedding(kernel_basis_31);  # becomes the 0 ideal

# Chart 2: u[2] = 1
kernel_basis_32 = ideal(D3, union(basis_blowup_3, [u[2] - 1]));
b_minimal_vars_32, b_minimal_ideal_32, b_elim_rules_32 = minimal_embedding(kernel_basis_32);
X_32 = spec(D3, b_minimal_ideal_32);
print(is_smooth(X_32));  # returns true

# Chart 3: u[3] = 1
kernel_basis_33 = ideal(D3, union(basis_blowup_3, [u[3] - 1]));
b_minimal_vars_33, b_minimal_ideal_33, b_elim_rules_33 = minimal_embedding(kernel_basis_33);
X_33 = spec(D3, b_minimal_ideal_33);
print(is_smooth(X_33));  # returns true

# Chart 4: u[4] = 1
kernel_basis_34 = ideal(D3, union(basis_blowup_3, [u[4] - 1]));
b_minimal_vars_34, b_minimal_ideal_34, b_elim_rules_34 = minimal_embedding(kernel_basis_34);
eval_b_min_34 = evaluate_variables(gens(b_minimal_ideal_34), [15]);
nonzero_eval_b_min_34 = nonzero_constant_indices(eval_b_min_34);  # Shows that there is a non-zero constant generator when u[1] is set to zero, proving that the first chart covers it
print(nonzero_eval_b_min_34);

# Chart 5: u[5] = 1
kernel_basis_35 = ideal(D3, union(basis_blowup_3, [u[5] - 1]));
b_minimal_vars_35, b_minimal_ideal_35, b_elim_rules_35 = minimal_embedding(kernel_basis_35);
X_35 = spec(D3, b_minimal_ideal_35);
print(is_smooth(X_35));  # returns true

# Chart 6: u[6] = 1
kernel_basis_36 = ideal(D3, union(basis_blowup_3, [u[6] - 1]));
b_minimal_vars_36, b_minimal_ideal_36, b_elim_rules_36 = minimal_embedding(kernel_basis_36);  # becomes the zero ideal

# Chart 7: u[7] = 1
kernel_basis_37 = ideal(D3, union(basis_blowup_3, [u[7] - 1]));
b_minimal_vars_37, b_minimal_ideal_37, b_elim_rules_37 = minimal_embedding(kernel_basis_37);
X_37 = spec(D3, b_minimal_ideal_37);
print(is_smooth(X_37));  # returns true

# Chart 8: u[8] = 1
kernel_basis_38 = ideal(D3, union(basis_blowup_3, [u[8] - 1]));
b_minimal_vars_38, b_minimal_ideal_38, b_elim_rules_38 = minimal_embedding(kernel_basis_38);
X_38 = spec(D3, b_minimal_ideal_38);
print(is_smooth(X_38));  # returns true

# Chart 9: u[9] = 1
kernel_basis_39 = ideal(D3, union(basis_blowup_3, [u[9] - 1]));
b_minimal_vars_39, b_minimal_ideal_39, b_elim_rules_39 = minimal_embedding(kernel_basis_39);
eval_b_min_39 = evaluate_variables(gens(b_minimal_ideal_39), [15]);
nonzero_eval_b_min_39 = nonzero_constant_indices(eval_b_min_39);  # Shows that there is a non-zero constant generator when u[1] is set to zero, proving that the first chart covers it
print(nonzero_eval_b_min_39);

# Chart 10: u[10] = 1
kernel_basis_310 = ideal(D3, union(basis_blowup_3, [u[10] - 1]));
b_minimal_vars_310, b_minimal_ideal_310, b_elim_rules_310 = minimal_embedding(kernel_basis_310);
eval_b_min_310 = evaluate_variables(gens(b_minimal_ideal_310), [19, 20, 22]);
nonzero_eval_b_min_310 = nonzero_constant_indices(eval_b_min_310);  # Shows that there is a non-zero constant generator when u[5], u[6] and u[8] are set to zero, proving that the union of these charts covers it
print(nonzero_eval_b_min_310);

# Chart 11: u[11] = 1
kernel_basis_311 = ideal(D3, union(basis_blowup_3, [u[11] - 1]));
b_minimal_vars_311, b_minimal_ideal_311, b_elim_rules_311 = minimal_embedding(kernel_basis_311);
eval_b_min_311 = evaluate_variables(gens(b_minimal_ideal_311), [15]);
nonzero_eval_b_min_311 = nonzero_constant_indices(eval_b_min_311);  # Shows that there is a non-zero constant generator when u[1] is set to zero, proving that the first chart covers it
print(nonzero_eval_b_min_311);

# Chart 12: u[12] = 1
kernel_basis_312 = ideal(D3, union(basis_blowup_3, [u[12] - 1]));
b_minimal_vars_312, b_minimal_ideal_312, b_elim_rules_312 = minimal_embedding(kernel_basis_312);
X_312 = spec(D3, b_minimal_ideal_312);
print(is_smooth(X_312));  # returns true

# Chart 13: u[13] = 1
kernel_basis_313 = ideal(D3, union(basis_blowup_3, [u[13] - 1]));
b_minimal_vars_313, b_minimal_ideal_313, b_elim_rules_313 = minimal_embedding(kernel_basis_313);
eval_b_min_313 = evaluate_variables(gens(b_minimal_ideal_313), [15, 22]);
nonzero_eval_b_min_313 = nonzero_constant_indices(eval_b_min_313);  # Shows that there is a non-zero constant generator when u[1] and u[8] are set to zero, proving that the union of these charts covers it
print(nonzero_eval_b_min_313);

# Chart 14: u[14] = 1
kernel_basis_314 = ideal(D3, union(basis_blowup_3, [u[14] - 1]));
b_minimal_vars_314, b_minimal_ideal_314, b_elim_rules_314 = minimal_embedding(kernel_basis_314);
eval_b_min_314 = evaluate_variables(gens(b_minimal_ideal_314), [22, 24, 27]);
nonzero_eval_b_min_314 = nonzero_constant_indices(eval_b_min_314);  # Shows that there is a non-zero constant generator when u[8], u[10] and u[13] are set to zero, proving that the union of these charts covers it
print(nonzero_eval_b_min_314);

# Chart 15: u[15] = 1
kernel_basis_315 = ideal(D3, union(basis_blowup_3, [u[15] - 1]));
b_minimal_vars_315, b_minimal_ideal_315, b_elim_rules_315 = minimal_embedding(kernel_basis_315);
eval_b_min_315 = evaluate_variables(gens(b_minimal_ideal_315), [19, 20, 22]);
nonzero_eval_b_min_315 = nonzero_constant_indices(eval_b_min_315);  # Shows that there is a non-zero constant generator when u[5], u[6] and u[8] are set to zero, proving that the union of these charts covers it
print(nonzero_eval_b_min_315);

# Chart 16: u[16] = 1
kernel_basis_316 = ideal(D3, union(basis_blowup_3, [u[16] - 1]));
b_minimal_vars_316, b_minimal_ideal_316, b_elim_rules_316 = minimal_embedding(kernel_basis_316);
eval_b_min_316 = evaluate_variables(gens(b_minimal_ideal_316), [31]);
nonzero_eval_b_min_316 = nonzero_constant_indices(eval_b_min_316);  # Shows that there is a non-zero constant generator when u[17]is to zero, proving that this chart covers it
print(nonzero_eval_b_min_316);

# Chart 17: u[17] = 1
kernel_basis_317 = ideal(D3, union(basis_blowup_3, [u[17] - 1]));
b_minimal_vars_317, b_minimal_ideal_317, b_elim_rules_317 = minimal_embedding(kernel_basis_317);
X_317 = spec(D3, b_minimal_ideal_317);
print(is_smooth(X_317));  # returns true


D6, z, b, y = polynomial_ring(M, :z => (1:8), :b => (1:6), :y => (1:4));
blowup_map_6 = hom(D6, Ct, blowup_vec_6);
relations_chart_6 = ideal(Ct, [embed_Ct(f) for f in gens(c_6_simplified_ideal)]);
kernel_blowup_map_6 = preimage(blowup_map_6, relations_chart_6);

# Chart 1: y[1] = 1
kernel_basis_61 = ideal(D6, union(gens(kernel_blowup_map_6), [y[1] - 1]));
b_minimal_vars_61, b_minimal_ideal_61, b_elim_rules_61 = minimal_embedding(kernel_basis_61);
X_61 = spec(D6, b_minimal_ideal_61);
print(is_smooth(X_61));  # returns true?

# Chart 2: y[2] = 1
kernel_basis_62 = ideal(D6, union(gens(kernel_blowup_map_6), [y[2] - 1]));
b_minimal_vars_62, b_minimal_ideal_62, b_elim_rules_62 = minimal_embedding(kernel_basis_62);
X_62 = spec(D6, b_minimal_ideal_62);
print(is_smooth(X_62));  # returns true?

# Chart 3: y[3] = 1
kernel_basis_63 = ideal(D6, union(gens(kernel_blowup_map_6), [y[3] - 1]));
b_minimal_vars_63, b_minimal_ideal_63, b_elim_rules_63 = minimal_embedding(kernel_basis_63);
X_63 = spec(D6, b_minimal_ideal_63);
print(is_smooth(X_63));  # returns true?

# Chart 4: y[4] = 1
kernel_basis_64 = ideal(D6, union(gens(kernel_blowup_map_6), [y[4] - 1]));
b_minimal_vars_64, b_minimal_ideal_64, b_elim_rules_64 = minimal_embedding(kernel_basis_64);
X_64 = spec(D6, b_minimal_ideal_64);
print(is_smooth(X_64));  # returns true?

# Computations of the fiber of the origin that were initially done to distinguish the two symplectic resolutions but they are no longer relevant
# Only the dimension of the fiber is relevant, but this is covered by the symmetric case in blowupsW2.jl
#=
# Finally, we compute the fiber of the blowup map over the origin
fiber_origin_basis_ideal = ideal(C, union(gens(basis_blowup_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
reduced_fiber_origin_basis_ideal = radical(fiber_origin_basis_ideal);

# The reduced fiber is used to compute the dimension of the fiber, which is used to show that the resolution is semi-small
X_fiber = spec(C, reduced_fiber_origin_basis_ideal);
components_fiber = irreducible_components(X_fiber);

# In order to compare the fibers of W_1 and W_2, we use the non-reduced fiber
X_nonred_fiber = spec(C, fiber_origin_basis_ideal);

min_basis_fiber = standard_basis(modulus(OO(X_nonred_fiber)), ordering=negdegrevlex(C));
X_min_basis_fiber = spec(C, ideal(C, min_basis_fiber));

sing_fiber = singular_locus(X_min_basis_fiber);
sing_ideal = modulus(OO(sing_fiber[1]));
reduced_sing_ideal = radical(sing_ideal);  # Shows that the reduced singular locus is not the whole fiber
=#