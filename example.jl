using Oscar
Qx, s = polynomial_ring(QQ, :s);

# For the rest of the computations, we only need the algebraic number sqrt(-3), hence we can pass to a smaller number field
M, q3 = number_field(s^2 + 3, "q3");

Bt, z, t = polynomial_ring(M, :z => (1:8), :t => (1:1));

rel_1 = q3*z[1]^3*z[5] - z[1]*z[3]*z[4] - 2*z[2]*z[6] - z[5]*z[8];
rel_2 = q3*z[1]^3*z[4] + z[1]*z[2]*z[5] - 2*z[3]*z[7] - z[4]*z[8];
rel_3 = z[1]*z[5]^2 + 2*z[4]*z[6] + z[3]*z[8];
rel_4 = z[1]*z[4]^2 - 2*z[5]*z[7] - z[2]*z[8];
rel_5 = -z[1]^4 + z[2]*z[3] - z[4]*z[5] - 3*q3*z[1]*z[8];
rel_6 = z[1]^2*z[4]*z[5] + q3*z[1]^3*z[8] + 4*z[6]*z[7] - z[8]^2;
rel_7 = q3*z[1]^2*z[3]*z[5] - 2*z[1]^3*z[6] - z[3]^2*z[4] + z[5]^3 - 6*q3*z[6]*z[8];
rel_8 = q3*z[1]^2*z[2]*z[4] - 2*z[1]^3*z[7] - z[4]^3 + z[2]^2*z[5] - 6*q3*z[7]*z[8];
rel_9 = 4*z[1]^2*z[4]*z[5] + q3*z[3]*z[4]^2 - q3*z[2]*z[5]^2 + 4*z[6]*z[7] + 8*z[8]^2;
relations = ideal(Bt, [rel_1, rel_2, rel_3, rel_4, rel_5, rel_6, rel_7, rel_8, rel_9]);

W_gen_1 = z[3]*z[7] + 2*z[4]*z[8];
W_gen_2 = z[2]*z[4] + 2*q3*z[1]*z[7];
W_gen_3 = z[2]*z[3] - 4*q3*z[1]*z[8];
W_gen_4 = z[2]^3 + 12*q3*z[7]^2;
W_gen_5 = z[1]*z[2]^2 - 6*z[4]*z[7];
W_gen_6 = z[1]^2*z[2] - q3*z[4]^2;
weildivisor = ideal(Bt, [W_gen_1, W_gen_2, W_gen_3, W_gen_4, W_gen_5, W_gen_6]);

# Define the vector used to construct the underlying map of the blowup along the Weil divisor W_2
blowup_vector = union([z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]] , gens(weildivisor*t[1]));

C, z, b = polynomial_ring(M, :z => (1:8), :b => (1:6));
blowup_map = hom(C, Bt, blowup_vector);
kernel_blowup_map = preimage(blowup_map, relations);
basis_blowup = standard_basis(kernel_blowup_map, ordering=negdegrevlex(C));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal = ideal(C, collect(basis_blowup));

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
                        break  # We found a variable to eliminate and thus can break the inner loop
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


# Chart 1: b[1] != 0
kernel_basis_1 = ideal(C, union(gens(basis_blowup_ideal), [b[1] - 1]));
b_minimal_vars_1, b_minimal_ideal_1, b_elim_rules_1 = minimal_embedding(kernel_basis_1);
F_1_b = free_resolution(b_minimal_ideal_1, algorithm=:mres);
