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

# Define the two Weil divisors that are blown up to establish the first partial resolution
W_gen_11 = z[2]*z[6] + 2*z[5]*z[8];
W_gen_12 = z[3]*z[5] + 2*q3*z[1]*z[6];
W_gen_13 = z[2]*z[3] - 4*q3*z[1]*z[8];
W_gen_14 = z[3]^3 - 12*q3*z[6]^2;
W_gen_15 = z[1]*z[3]^2 + 6*z[5]*z[6];
W_gen_16 = z[1]^2*z[3] + q3*z[5]^2;
weildivisor_1 = ideal(Bt, [W_gen_11, W_gen_12, W_gen_13, W_gen_14, W_gen_15, W_gen_16]);

W_gen_21 = z[3]*z[7] + 2*z[4]*z[8];
W_gen_22 = z[2]*z[4] + 2*q3*z[1]*z[7];
W_gen_23 = z[2]*z[3] - 4*q3*z[1]*z[8];
W_gen_24 = z[2]^3 + 12*q3*z[7]^2;
W_gen_25 = z[1]*z[2]^2 - 6*z[4]*z[7];
W_gen_26 = z[1]^2*z[2] - q3*z[4]^2;
weildivisor_2 = ideal(Bt, [W_gen_21, W_gen_22, W_gen_23, W_gen_24, W_gen_25, W_gen_26]);

# Define the vector used to construct the underlying map of the blowup along the Weil divisor W_1
blowup_vector_1 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]] , gens(weildivisor_1*t[1]));

# Define the vector used to construct the underlying map of the blowup along the Weil divisor W_1
blowup_vector_2 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]] , gens(weildivisor_2*t[1]));

# Define the ring C_1=B[b_1,...,b_6] used to compute the rees algebra of the blowups along the Weil divisor W_1
C_1, z, b = polynomial_ring(M, :z => (1:8), :b => (1:6));
blowup_map_1 = hom(C_1, Bt, blowup_vector_1);
kernel_blowup_map_1 = preimage(blowup_map_1, relations);
basis_blowup_1 = standard_basis(kernel_blowup_map_1, ordering=negdegrevlex(C_1));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_1 = ideal(C_1, collect(basis_blowup_1));

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

# In the following, we only consider charts in detail if they are involved in the second blowup

# Chart 1: b[1] = 1
chart_11_ideal = ideal(C_1, union(basis_blowup_1, [b[1]-1]));
#b_minimal_vars_11, b_minimal_ideal_11, b_elim_rules_11 = minimal_embedding(chart_11_ideal);


# Chart 2: b[2] = 1
chart_12_ideal = ideal(C_1, union(basis_blowup_1, [b[2]-1]));
b_minimal_vars_12, b_minimal_ideal_12, b_elim_rules_12 = minimal_embedding(chart_12_ideal);

generator_12 = 12*z[1]^2 + 3*q3*b[3]*b[4] - 4*q3*z[1]*b[5]*b[6] - b[5]^2*b[6]^2 + b[4]*b[6]^3;  # Simplified generator found by inspecting the standard basis
c_12_simplified_ideal = ideal(C_1, generator_12);
print(b_minimal_ideal_12 == c_12_simplified_ideal);  # Shows that the ideal can be simplified to only one generator

X_12 = spec(C_1, c_12_simplified_ideal);
sing_12 = singular_locus(X_12);
sing_12_ideal = modulus(OO(sing_12[1]));

chart_12_reconstructed_ideal = ideal(C_1, union(gens(b_minimal_ideal_12), [k - b_elim_rules_12[k] for k in keys(b_elim_rules_12)]));
print(chart_12_ideal == chart_12_reconstructed_ideal);  # Shows that after adding back in the elimination rules, we get the same variety, the ideal defining the singular locus remains untouched by the elimination rules


# Chart 3: b[3] = 1
chart_13_ideal = ideal(C_1, union(basis_blowup_1, [b[3]-1]));
b_minimal_vars_13, b_minimal_ideal_13, b_elim_rules_13 = minimal_embedding(chart_13_ideal);
generator_13_1 = gens(b_minimal_ideal_13)[11];
generator_13_2 = gens(b_minimal_ideal_13)[16];
generator_13_3 = gens(b_minimal_ideal_13)[19];
generator_13_4 = gens(b_minimal_ideal_13)[20];
generator_13_5 = gens(b_minimal_ideal_13)[23];
c_13_simplified_ideal = ideal(C_1, [generator_13_1, generator_13_2, generator_13_3, generator_13_4, generator_13_5]);
print(b_minimal_ideal_13 == c_13_simplified_ideal);  # Shows that the ideal can be simplified to only five generators

X_13 = spec(C_1, c_13_simplified_ideal);
sing_13 = singular_locus(X_13);
sing_13_ideal = modulus(OO(sing_13[1]));

chart_13_reconstructed_ideal = ideal(C_1, union(gens(b_minimal_ideal_13), [k - b_elim_rules_13[k] for k in keys(b_elim_rules_13)]));
print(chart_13_ideal == chart_13_reconstructed_ideal);  # Shows that after adding back in the elimination rules, we get the same variety, the ideal defining the singular locus remains untouched by the elimination rules


# Chart 4: b[4] = 1
chart_14_ideal = ideal(C_1, union(basis_blowup_1, [b[4]-1]));
#b_minimal_vars_14, b_minimal_ideal_14, b_elim_rules_14 = minimal_embedding(kernel_basis_14);


# Chart 5: b[5] = 1
chart_15_ideal = ideal(C_1, union(basis_blowup_1, [b[5]-1]));
#b_minimal_vars_15, b_minimal_ideal_15, b_elim_rules_15 = minimal_embedding(kernel_basis_15);


# Chart 6: b[6] = 1
chart_16_ideal = ideal(C_1, union(basis_blowup_1, [b[6]-1]));
b_minimal_vars_16, b_minimal_ideal_16, b_elim_rules_16 = minimal_embedding(chart_16_ideal);
gen_16_1 = z[5] + 2*q3*b[1]*b[5] + 4*z[4]*b[2]^2 - q3*z[5]*b[2]^2*b[3] - 2*b[1]*b[2]^2*b[3]*b[5] - b[2]*b[3]^2*b[5]^2;
gen_16_2 =  3*z[4] - 2*q3*z[5]*b[3] + 12*b[1]^2*b[2] + 4*q3*z[4]*b[2]^2*b[3] + 3*z[5]*b[2]^2*b[3]^2 - 2*q3*b[1]*b[2]^2*b[3]^2*b[5] - q3*b[2]*b[3]^3*b[5]^2;
c_16_simplified_ideal = ideal(C_1, [gen_16_1, gen_16_2]);  # Simplified generators found by inspecting the standard basis
print(b_minimal_ideal_16 == c_16_simplified_ideal);  # Shows that the ideal can be simplified to only two generators

X_16 = spec(C_1, c_16_simplified_ideal);
sing_16 = singular_locus(X_16);
sing_16_ideal = modulus(OO(sing_16[1]));

chart_16_reconstructed_ideal = ideal(C_1, union(gens(b_minimal_ideal_16), [k - b_elim_rules_16[k] for k in keys(b_elim_rules_16)]));
print(chart_16_ideal == chart_16_reconstructed_ideal);  # Shows that after adding back in the elimination rules, we get the same variety, the ideal defining the singular locus remains untouched by the elimination rules

# Compute the second blowup in the charts admitting singularities to get the resolution
min_sing_12 = standard_basis(sing_12_ideal, ordering=negdegrevlex(C_1));
min_sing_13 = standard_basis(sing_13_ideal, ordering=negdegrevlex(C_1));
min_sing_16 = gens(radical(sing_16_ideal));

C_1t, z, b, t = polynomial_ring(M, :z => (1:8), :b => (1:6), :t => (1:1));
embed_C_1t = hom(C_1, C_1t, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]]);

blowup_vec_12 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_C_1t(f)*t[1] for f in min_sing_12]);
blowup_vec_13 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_C_1t(f)*t[1] for f in min_sing_13]);
blowup_vec_16 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6]], [embed_C_1t(f)*t[1] for f in min_sing_16]);

# We compute the blowups locally above the covering charts having a singularity by using the Rees algebra as before
D_12, z, b, w = polynomial_ring(M, :z => (1:8), :b => (1:6), :w => (1:3));
blowup_map_12 = hom(D_12, C_1t, blowup_vec_12);
relations_chart_12 = ideal(C_1t, [embed_C_1t(f) for f in gens(chart_12_reconstructed_ideal)]);
kernel_blowup_map_12 = preimage(blowup_map_12, relations_chart_12);

chart_121_ideal = ideal(D_12, union(gens(kernel_blowup_map_12), [w[1]-1]));
zero_fiber_121 = ideal(D_12, union(gens(chart_121_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_121));
chart_122_ideal = ideal(D_12, union(gens(kernel_blowup_map_12), [w[2]-1]));
zero_fiber_122 = ideal(D_12, union(gens(chart_122_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_122));
chart_123_ideal = ideal(D_12, union(gens(kernel_blowup_map_12), [w[3]-1]));
zero_fiber_123 = ideal(D_12, union(gens(chart_123_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_123));

D_13, z, b, u = polynomial_ring(M, :z => (1:8), :b => (1:6), :u => (1:17));
blowup_map_13 = hom(D_13, C_1t, blowup_vec_13);
relations_chart_13 = ideal(C_1t, [embed_C_1t(f) for f in gens(chart_13_reconstructed_ideal)]);
#kernel_blowup_map_13 = preimage(blowup_map_13, relations_chart_13);  # inefficient

D_16, z, b, y = polynomial_ring(M, :z => (1:8), :b => (1:6), :y => (1:4));
blowup_map_16 = hom(D_16, C_1t, blowup_vec_16);
relations_chart_16 = ideal(C_1t, [embed_C_1t(f) for f in gens(chart_16_reconstructed_ideal)]);
#kernel_blowup_map_16 = preimage(blowup_map_16, relations_chart_16);  # inefficient

# Define the ring C=B[c_1,...,c_6] used to compute the rees algebra of the blowups along the Weil divisor W_2
C_2, z, c = polynomial_ring(M, :z => (1:8), :c => (1:6));
blowup_map_2 = hom(C_2, Bt, blowup_vector_2);
kernel_blowup_map_2 = preimage(blowup_map_2, relations);
basis_blowup_2 = standard_basis(kernel_blowup_map_2, ordering=negdegrevlex(C_2));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_2 = ideal(C_2, collect(basis_blowup_2));

# Chart 1: c[1] = 1
chart_21_ideal = ideal(C_2, union(basis_blowup_2, [c[1]-1]));
#b_minimal_vars_21, b_minimal_ideal_21, b_elim_rules_21 = minimal_embedding(chart_21_ideal);


# Chart 2: c[2] = 1
chart_22_ideal = ideal(C_2, union(basis_blowup_2, [c[2]-1]));
b_minimal_vars_22, b_minimal_ideal_22, b_elim_rules_22 = minimal_embedding(chart_22_ideal);

generator_22 = 12*z[1]^2 + 4*q3*z[1]*c[5]*c[6] - 3*q3*c[3]*c[4] + c[4]*c[6]^3 - c[5]^2*c[6]^2;  # Simplified generator found by inspecting the standard basis
c_22_simplified_ideal = ideal(C_2, generator_22);
print(b_minimal_ideal_22 == c_22_simplified_ideal)

X_22 = spec(C_2, c_22_simplified_ideal);
sing_22 = singular_locus(X_22);
sing_22_ideal = modulus(OO(sing_22[1]));

chart_22_reconstructed_ideal = ideal(C_2,union(gens(b_minimal_ideal_22), [k - b_elim_rules_22[k] for k in keys(b_elim_rules_22)]));
print(chart_22_ideal == chart_22_reconstructed_ideal);


# Chart 3: c[3] = 1
chart_23_ideal = ideal(C_2, union(basis_blowup_2, [c[3]-1]));
b_minimal_vars_23, b_minimal_ideal_23, b_elim_rules_23 = minimal_embedding(chart_23_ideal);
generator_23_1 = gens(b_minimal_ideal_23)[12];
generator_23_2 = gens(b_minimal_ideal_23)[15];
generator_23_3 = gens(b_minimal_ideal_23)[19];
generator_23_4 = gens(b_minimal_ideal_23)[21];
generator_23_5 = gens(b_minimal_ideal_23)[22];
c_23_simplified_ideal = ideal(C_2, [generator_23_1, generator_23_2, generator_23_3, generator_23_4, generator_23_5]);
print(b_minimal_ideal_23 == c_23_simplified_ideal);

X_23 = spec(C_2, c_23_simplified_ideal);
sing_23 = singular_locus(X_23);
sing_23_ideal = modulus(OO(sing_23[1]));

chart_23_reconstructed_ideal = ideal(C_2, union(gens(b_minimal_ideal_23), [k - b_elim_rules_23[k] for k in keys(b_elim_rules_23)]));
print(chart_23_ideal == chart_23_reconstructed_ideal);


# Chart 4: c[4] = 1
chart_24_ideal = ideal(C_2, union(basis_blowup_2, [c[4]-1]));
#b_minimal_vars_24, b_minimal_ideal_24, b_elim_rules_24 = minimal_embedding(chart_24_ideal);


# Chart 5: c[5] = 1
chart_25_ideal = ideal(C_2, union(basis_blowup_2, [c[5]-1]));
#b_minimal_vars_25, b_minimal_ideal_25, b_elim_rules_25 = minimal_embedding(chart_25_ideal);


# Chart 6: c[6] = 1
chart_26_ideal = ideal(C_2, union(basis_blowup_2, [c[6]-1]));
b_minimal_vars_26, b_minimal_ideal_26, b_elim_rules_26 = minimal_embedding(chart_26_ideal);

gen_6_1 = z[4] + 2*q3*c[1]*c[5] + 4*z[5]*c[2]^2 + q3*z[4]*c[2]^2*c[3] + 2*c[1]*c[2]^2*c[3]*c[5] - c[2]*c[3]^2*c[5]^2;
gen_6_2 =  3*z[5] + 2*q3*z[4]*c[3] + 12*c[1]^2*c[2] - 4*q3*z[5]*c[2]^2*c[3] + 3*z[4]*c[2]^2*c[3]^2 - 2*q3*c[1]*c[2]^2*c[3]^2*c[5] + q3*c[2]*c[3]^3*c[5]^2;
c_26_simplified_ideal = ideal(C_2, [gen_6_1, gen_6_2]);  # Simplified generators found by inspecting the standard basis
print(b_minimal_ideal_26 == c_26_simplified_ideal);

X_26 = spec(C_2, c_26_simplified_ideal);
sing_26 = singular_locus(X_26);
sing_26_ideal = modulus(OO(sing_26[1]));

chart_26_reconstructed_ideal = ideal(C_2, union(gens(b_minimal_ideal_26), [k - b_elim_rules_26[k] for k in keys(b_elim_rules_26)]));
print(chart_26_ideal == chart_26_reconstructed_ideal); 

# Compute the second blowup in the charts admitting singularities to get the resolution
min_sing_22 = standard_basis(sing_22_ideal, ordering=negdegrevlex(C_2));
min_sing_23 = standard_basis(sing_23_ideal, ordering=negdegrevlex(C_2));
min_sing_26 = gens(radical(sing_26_ideal));

C_2t, z, c, t = polynomial_ring(M, :z => (1:8), :c => (1:6), :t => (1:1));
embed_C_2t = hom(C_2, C_2t, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], c[1], c[2], c[3], c[4], c[5], c[6]]);

blowup_vec_22 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], c[1], c[2], c[3], c[4], c[5], c[6]], [embed_C_2t(f)*t[1] for f in min_sing_22]);
blowup_vec_23 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], c[1], c[2], c[3], c[4], c[5], c[6]], [embed_C_2t(f)*t[1] for f in min_sing_23]);
blowup_vec_26 = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], c[1], c[2], c[3], c[4], c[5], c[6]], [embed_C_2t(f)*t[1] for f in min_sing_26]);

# We compute the blowups locally above the covering charts having a singularity by using the Rees algebra as before
D_22, z, c, w = polynomial_ring(M, :z => (1:8), :c => (1:6), :w => (1:3));
blowup_map_22 = hom(D_22, C_2t, blowup_vec_22);
relations_chart_22 = ideal(C_2t, [embed_C_2t(f) for f in gens(chart_22_reconstructed_ideal)]);
kernel_blowup_map_22 = preimage(blowup_map_22, relations_chart_22);

chart_221_ideal = ideal(D_22, union(gens(kernel_blowup_map_22), [w[1]-1]));
zero_fiber_221 = ideal(D_22, union(gens(chart_221_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_221));
chart_222_ideal = ideal(D_22, union(gens(kernel_blowup_map_22), [w[2]-1]));
zero_fiber_222 = ideal(D_22, union(gens(chart_222_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_222));
chart_223_ideal = ideal(D_22, union(gens(kernel_blowup_map_22), [w[3]-1]));
zero_fiber_223 = ideal(D_22, union(gens(chart_223_ideal), [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
print(dim(zero_fiber_223));

D_23, z, c, u = polynomial_ring(M, :z => (1:8), :c => (1:6), :u => (1:17));
blowup_map_23 = hom(D_23, C_2t, blowup_vec_23);
relations_chart_23 = ideal(C_2t, [embed_C_2t(f) for f in gens(chart_23_reconstructed_ideal)]);
#kernel_blowup_map_23 = preimage(blowup_map_23, relations_chart_23);  # inefficient

D_26, z, c, y = polynomial_ring(M, :z => (1:8), :c => (1:6), :y => (1:4));
blowup_map_26 = hom(D_26, C_2t, blowup_vec_26);
relations_chart_26 = ideal(C_2t, [embed_C_2t(f) for f in gens(chart_26_reconstructed_ideal)]);
#kernel_blowup_map_26 = preimage(blowup_map_26, relations_chart_26);

# Fiber products of the affine charts corresponding to the blowups of the second charts
P, z, b, c, = polynomial_ring(M, :z => (1:8), :b => (1:9), :c => (1:9));
embed_D12 = hom(D_12, P, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9]]);
embed_D22 = hom(D_22, P, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9]]);

prod_1 = ideal(P, union([embed_D12(f) for f in gens(chart_121_ideal)], [embed_D22(f) for f in gens(chart_221_ideal)]));
prod_2 = ideal(P, union([embed_D12(f) for f in gens(chart_122_ideal)], [embed_D22(f) for f in gens(chart_222_ideal)]));
prod_3 = ideal(P, union([embed_D12(f) for f in gens(chart_123_ideal)], [embed_D22(f) for f in gens(chart_223_ideal)]));

# Consider the fiber of the origin in these fiber products
prod_fiber_1 = ideal(P, union([embed_D12(f) for f in gens(chart_121_ideal)], [embed_D22(f) for f in gens(chart_221_ideal)], [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
prod_fiber_2 = ideal(P, union([embed_D12(f) for f in gens(chart_122_ideal)], [embed_D22(f) for f in gens(chart_222_ideal)], [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));
prod_fiber_3 = ideal(P, union([embed_D12(f) for f in gens(chart_123_ideal)], [embed_D22(f) for f in gens(chart_223_ideal)], [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]));

X_prod_1 = spec(P, prod_1);
X_prod_2 = spec(P, prod_2);
X_prod_3 = spec(P, prod_3);

# The fiber of the graph has to consist of some of the irreducible components of the fiber product
X_prod_fiber_1 = spec(P, prod_fiber_1);
comp_fiber_1 = irreducible_components(X_prod_fiber_1);
print([dim(c) for c in comp_fiber_1]);  # all of dimension 4

X_prod_fiber_2 = spec(P, prod_fiber_2);
comp_fiber_2 = irreducible_components(X_prod_fiber_2);
print([dim(c) for c in comp_fiber_2]);  # all of dimension 4

X_prod_fiber_3 = spec(P, prod_fiber_3);
comp_fiber_3 = irreducible_components(X_prod_fiber_3);
print([dim(c) for c in comp_fiber_3]);  # all of dimension 4

origin_ideal = ideal(P, [z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]]);

# The closure of the graph is computed via saturation, 
# Note that none of the following computations ran through in a reasonable running time whence the results are unclear
graph_closure_1 = saturation(prod_1, origin_ideal);
graph_closure_2 = saturation(prod_2, origin_ideal);
graph_closure_3 = saturation(prod_3, origin_ideal);

# Compute the central fiber of the graph closure
graph_fiber_1 = graph_closure_1 + origin_ideal;
graph_fiber_2 = graph_closure_2 + origin_ideal;
graph_fiber_3 = graph_closure_3 + origin_ideal;

# Check the dimensions of the true graph fibers
print("Dimension of Graph Fiber 1: ", dim(graph_fiber_1), "\n");
print("Dimension of Graph Fiber 2: ", dim(graph_fiber_2), "\n");
print("Dimension of Graph Fiber 3: ", dim(graph_fiber_3), "\n");