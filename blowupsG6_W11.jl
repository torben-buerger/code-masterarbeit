using Oscar
Qx, s = polynomial_ring(QQ, :s);

# For the rest of the computations, we only need the algebraic number sqrt(-3), hence we can pass to a smaller number field
M, q3 = number_field(s^2 + 3, "q3");

# In order to compute the blowups along the previously computed varieties belonging to delta_11, delta_21, and delta_11_21, we first define a polynomial ring that represents the coordinate ring of the quotient variety
# We add the variable t to compute the Rees algebra of the ideal of the variety we want to blow up
Bt, y, t = polynomial_ring(M, :y => (1:10), :t => (1:1));

# The coordinate ring of the quotient variety corresponds to the quotient of Bt by the relation ideal
relations_ideal = ideal(Bt, [-3*y[1]^6 - 10*q3*y[1]^3*y[4] + 3*y[1]^2*y[2]*y[3] + q3*y[2]*y[6] + 3*q3*y[3]*y[5] + 9*y[4]^2,
  y[1]^5*y[2] - 3*q3*y[1]^3*y[5] + 3*q3*y[1]^2*y[2]*y[4] - y[1]*y[2]^2*y[3] + 2*q3*y[3]*y[7] + 3*y[4]*y[5],
  y[1]^4*y[2]^2 + 2*q3*y[1]^3*y[7] - 3*q3*y[1]^2*y[2]*y[5] + 3*q3*y[1]*y[2]^2*y[4] - y[2]^3*y[3] - 18*y[4]*y[7] + 9*y[5]^2,
  y[1]^5*y[3] - q3*y[1]^3*y[6] + 3*q3*y[1]^2*y[3]*y[4] - y[1]*y[2]*y[3]^2 + 2*y[2]*y[8] + y[4]*y[6],
  8*q3*y[1]^5*y[4] - 3*y[1]^4*y[2]*y[3] + q3*y[1]^2*y[2]*y[6] + 3*q3*y[1]^2*y[3]*y[5] - 72*y[1]^2*y[4]^2 - 18*q3*y[1]*y[2]*y[3]*y[4] + 3*y[2]^2*y[3]^2 + 9*y[5]*y[6],
  y[1]^4*y[3]^2 + 2*y[1]^3*y[8] - q3*y[1]^2*y[3]*y[6] + 3*q3*y[1]*y[3]^2*y[4] - y[2]*y[3]^3 + 6*q3*y[4]*y[8] + y[6]^2,
  -18*q3*y[1]^3*y[7] + 9*q3*y[1]^2*y[2]*y[5] - 3*q3*y[1]*y[2]^2*y[4] + 2*y[2]^3*y[3] - 2*y[3]*y[10] + 18*y[4]*y[7],
  2*y[1]^3*y[2]^3 - 2*y[1]^3*y[10] - 18*q3*y[1]^2*y[2]*y[7] + 9*q3*y[1]*y[2]^2*y[5] + 3*q3*y[2]^3*y[4] - 6*q3*y[4]*y[10] + 54*y[5]*y[7],
  12*q3*y[1]^5*y[5] - 6*q3*y[1]^2*y[3]*y[7] - 36*y[1]^2*y[4]*y[5] - q3*y[1]*y[2]^2*y[6] - 6*q3*y[1]*y[2]*y[3]*y[5] + q3*y[2]^2*y[3]*y[4] + 6*y[6]*y[7],
  q3*y[2]^3*y[5] - q3*y[5]*y[10] + 18*y[7]^2,
  -18*q3*y[1]^3*y[8] - 9*y[1]^2*y[3]*y[6] + 9*y[1]*y[3]^2*y[4] + 2*q3*y[2]*y[3]^3 - 2*q3*y[2]*y[9] + 18*y[4]*y[8],
  -4*y[1]^5*y[6] - 2*q3*y[1]^2*y[2]*y[8] - 4*q3*y[1]^2*y[4]*y[6] + 2*y[1]*y[2]*y[3]*y[6] + 3*y[1]*y[3]^2*y[5] - y[2]*y[3]^2*y[4] + 6*y[5]*y[8],
  2*q3*y[1]^3*y[3]^3 - 2*q3*y[1]^3*y[9] - 18*q3*y[1]^2*y[3]*y[8] - 9*y[1]*y[3]^2*y[6] - 9*y[3]^3*y[4] + 18*y[4]*y[9] + 18*y[6]*y[8],
  -3*q3*y[1]^2*y[5]*y[6] + q3*y[1]*y[2]*y[4]*y[6] + 3*q3*y[1]*y[3]*y[4]*y[5] - q3*y[2]*y[3]*y[4]^2 + 12*y[7]*y[8],
  - q3*y[3]^3*y[6] + q3*y[6]*y[9] + 18*y[8]^2,
  q3*y[1]^2*y[6]^2 - 2*q3*y[1]*y[3]*y[4]*y[6] - 2*y[3]^3*y[5] + q3*y[3]^2*y[4]^2 + 2*y[5]*y[9],
  72*y[1]^2*y[5]*y[8] - 3*q3*y[1]*y[2]*y[6]^2 - 9*q3*y[1]*y[3]*y[5]*y[6] - 24*y[1]*y[4]^2*y[6] + 3*q3*y[2]*y[3]*y[4]*y[6] - 4*y[3]^3*y[7] + 9*q3*y[3]^2*y[4]*y[5] + 24*y[3]*y[4]^3 + 4*y[7]*y[9],
  27*q3*y[1]^2*y[5]^2 - 18*q3*y[1]*y[2]*y[4]*y[5] - 2*y[2]^3*y[6] + 3*q3*y[2]^2*y[4]^2 + 2*y[6]*y[10],
  -72*y[1]^2*y[6]*y[7] + 27*y[1]*y[2]*y[5]*y[6] + 81*y[1]*y[3]*y[5]^2 - 72*q3*y[1]*y[4]^2*y[5] - 4*y[2]^3*y[8] - 9*y[2]^2*y[4]*y[6] - 27*y[2]*y[3]*y[4]*y[5] + 24*q3*y[2]*y[4]^3 + 4*y[8]*y[10], 
  -3290101285742342628191948215568900615539620384719577*y[1]^4*y[5]*y[6] + 30232015654944220871070857400173055270528000*y[1]^3*y[2]^2*y[8] - 367565213137716042200644235006599469144736404400429*y[1]^3*y[2]*y[4]*y[6]
  + 6046403130988844174214171480034611054105600*q3*y[1]^3*y[3]^2*y[7] - 6046403130988844174214171480034611054105600*q3*y[1]^2*y[2]^2*y[3]*y[6] + 435341025431196780543420346562491995895603200*q3*y[1]^2*y[7]*y[8] 
  - 1961765378255416525815782309542616299237080243702984*q3*y[1]*y[2]*y[5]*y[8] + 228637438148936169846427640834555114170169574692820*y[1]*y[3]*y[6]*y[7] - 1096700210910268160465592466812793590600542180438259*q3*y[1]*y[4]*y[5]*y[6]
  - 219306206037531022865697756718202521191546433397360*q3*y[2]^2*y[4]*y[8] + 38605887453810897682776120840324117262236355513989*y[2]^2*y[6]^2 + 1098199285706500867598067798074718076054956886430616*y[2]*y[3]*y[5]*y[6]
  - 225470863634249083049895071860360162944911912456647*q3*y[2]*y[4]^2*y[6] - 2015467710329614724738057160011537018035200*y[3]^3*y[10] - 654921048974828229261183332191220675988471082344966*y[3]^2*y[4]*y[7]
  + 2947144870035204523649218306661237172804743459665947*y[3]^2*y[5]^2 - 2619684631240338348241513872185229266445880224983064*q3*y[3]*y[4]^2*y[5]
  - 870682050862393561086840693124983991791206400*y[4]^4 + 2015467710329614724738057160011537018035200*y[9]*y[10]]);

# Define the genarators of the ideal W_delta_11 in which we want to compute the blowup via the Rees algebra
w_11_1 = y[7];
w_11_2 = -y[2]^3 + y[10];
w_11_3 = -3*y[1]*y[5] + y[2]*y[4];
w_11_4 = 4*y[1]^3*y[4] - y[2]*y[6] - 3*y[3]*y[5] + 4*q3*y[4]^2;
# Define the vector used to construct the underlying map of the blowup along the divisor W_delta_11
blowup_vector_11 = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], t[1]*w_11_1, t[1]*w_11_2, t[1]*w_11_3, t[1]*w_11_4];

# Define the ring C=B[b_1,...,b_4] used to compute the rees algebra of the blowups in W_delta_11 and W_delta_21, which both can be embedded into this relative projective space as the ideals along which we blow up are generated by 4 elements
C, y, b = polynomial_ring(M, :y => (1:10), :b => (1:4));
Aux_large_C, y, t, b = polynomial_ring(M, :y => (1:10), :t => (1:1), :b => (1:4));
embed_B_1 = hom(Bt, Aux_large_C, union([y[i] for i in 1:10], [t[1]]));
embed_C = hom(C, Aux_large_C, union([y[i] for i in 1:10], [b[i] for i in 1:4]));

# Compute the blowup along W_delta_11
blowup_map_11 = hom(C, Bt, blowup_vector_11);
kernel_blowup_map_11 = preimage(blowup_map_11, relations_ideal);
basis_blowup_11 = standard_basis(kernel_blowup_map_11, ordering=negdegrevlex(C));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_11 = ideal(C, collect(basis_blowup_11));
# Elimination theory approach to get the ideal of the blowup variety in the relative projective space P^4_C
blowup_ideal_11 = ideal(Aux_large_C, union([embed_B_1(rel) for rel in gens(relations_ideal)], [b[i] - embed_B_1(blowup_vector_11[10+i])*t[1] for i in 1:4]));
Rees_rel_11 = eliminate(blowup_ideal_11, [t[1]]);  # Get a presentation of the Rees algebra which corresponds to the blowup in delta_11
Rees_basis_11 = standard_basis(Rees_rel_11, ordering=negdegrevlex(Aux_large_C));  # Get a smaller generating set of the blowup ideal
Rees_basis_ideal_11 = ideal(Aux_large_C, Rees_basis_11);
print(embed_C(basis_blowup_ideal_11) == Rees_basis_ideal_11);  # Returns true


# Check that the blowup ideal is homogeneous with respect to the b-variables, i.e., in the ring B[b_1,...,b_4], so the y-variables have degree 0
# Note that this is done by evaluating all z-variables at 1 and checking homogeneity of the resulting polynomials in the b-variables because otherwise we would need to work 
# in a graded polynomial ring which does not allow the non-homogeneous generators needed to analyze the charts later on
b_part_kernel = [evaluate(gens(basis_blowup_11)[i], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) for i in eachindex(gens(basis_blowup_ideal_11))];
check_homogeneity = [is_homogeneous(b_part_kernel[i]) for i in eachindex(b_part_kernel)];  # Check that all generators are homogeneous with respect to the b-variables
print(all(check_homogeneity));

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


C, y, b = polynomial_ring(M, :y => (1:10), :b => (1:4));

# In the following we will investigate the specific charts based on their ideals. We will take both generating the, the original large one and the one given by the standard_basis
# Chart 1: b[1] = 1
kernel_1 = ideal(C, union(gens(basis_blowup_ideal_11), [b[1] - 1]));
minimal_vars_1, minimal_ideal_1, elim_rules_1 = minimal_embedding(kernel_1);
print(minimal_vars_1);
large_kernel_1 = ideal(C, union(gens(kernel_blowup_map_11), [b[1] - 1]));
large_minimal_vars_1, large_minimal_ideal_1, large_elim_rules_1 = minimal_embedding(large_kernel_1);

# Chart 2: b[2] = 1
kernel_2 = ideal(C, union(gens(basis_blowup_ideal_11), [b[2] - 1]));
minimal_vars_2, minimal_ideal_2, elim_rules_2 = minimal_embedding(kernel_2);
print(minimal_vars_2);
large_kernel_2 = ideal(C, union(gens(kernel_blowup_map_11), [b[2] - 1]));
large_minimal_vars_2, large_minimal_ideal_2, large_elim_rules_2 = minimal_embedding(large_kernel_2);

# Chart 3: b[3] = 1
kernel_3 = ideal(C, union(gens(basis_blowup_ideal_11), [b[3] - 1]));
minimal_vars_3, minimal_ideal_3, elim_rules_3 = minimal_embedding(kernel_3);
print(minimal_vars_1);
large_kernel_3 = ideal(C, union(gens(kernel_blowup_map_11), [b[3] - 1]));
large_minimal_vars_3, large_minimal_ideal_3, large_elim_rules_3 = minimal_embedding(large_kernel_3);

# Chart 1: b[4] = 1
kernel_4 = ideal(C, union(gens(basis_blowup_ideal_11), [b[4] - 1]));
minimal_vars_4, minimal_ideal_4, elim_rules_4 = minimal_embedding(kernel_4);
print(minimal_vars_4);
large_kernel_4 = ideal(C, union(gens(kernel_blowup_map_11), [b[4] - 1]));
large_minimal_vars_4, large_minimal_ideal_4, large_elim_rules_4 = minimal_embedding(large_kernel_4);

min_basis_1 = standard_basis(minimal_ideal_1, ordering=negdegrevlex(C));
print(length(min_basis_1));
b_minimal_ideal_1 = ideal(C, min_basis_1);
min_basis_2 = standard_basis(minimal_ideal_2, ordering=negdegrevlex(C));
print(length(min_basis_2));
b_minimal_ideal_2 = ideal(C, min_basis_2);
min_basis_3 = standard_basis(minimal_ideal_3, ordering=negdegrevlex(C));
print(length(min_basis_3));
b_minimal_ideal_3 = ideal(C, min_basis_3);
min_basis_4 = standard_basis(minimal_ideal_4, ordering=negdegrevlex(C));
print(length(min_basis_4));
b_minimal_ideal_4 = ideal(C, min_basis_4);


# Use the same functions as for G_4 to check, whether some chart is contained in another one
function evaluate_variables(gen_list, idx)  # Note that the indices of b[1],...,b[4] are 11,...,14 in C
    zeros_vec = zeros(Int, length(idx))
    return [evaluate(gen_list[i], idx, zeros_vec) for i in eachindex(gen_list)]
end

function nonzero_constant_indices(gen_list)  # Returns the indices of non-zero constant generators, used to detect constant generators after evaluation
    return [i for i in eachindex(gen_list)
            if is_constant(gen_list[i]) && gen_list[i] != 0]
end

eval_1 = evaluate_variables(gens(minimal_ideal_1), [12, 13, 14]);
nonzero_eval_gen_1 = nonzero_constant_indices(eval_1);  # Shows that there is no non-zero constant generator when all of b[2], b[3], and b[4] are set to 0
print(nonzero_eval_gen_1);
large_eval_1 = evaluate_variables(gens(large_minimal_ideal_1), [12, 13, 14]);
nonzero_large_eval_gen_1 = nonzero_constant_indices(large_eval_1);  
print(nonzero_large_eval_gen_1);

eval_2 = evaluate_variables(gens(minimal_ideal_2), [11, 13, 14]);
nonzero_eval_gen_2 = nonzero_constant_indices(eval_2);  # Shows that there is no non-zero constant generator when all of b[1], b[3], and b[4] are set to 0
print(nonzero_eval_gen_2);
large_eval_2 = evaluate_variables(gens(large_minimal_ideal_2), [11, 13, 14]);
nonzero_large_eval_gen_2 = nonzero_constant_indices(large_eval_2);  
print(nonzero_large_eval_gen_2);

eval_3 = evaluate_variables(gens(minimal_ideal_3), [11, 12, 14]);
nonzero_eval_gen_3 = nonzero_constant_indices(eval_3);  # Shows that there is no non-zero constant generator when all of b[1], b[2], and b[4] are set to 0
print(nonzero_eval_gen_3);
large_eval_3 = evaluate_variables(gens(large_minimal_ideal_3), [11, 12, 14]);
nonzero_large_eval_gen_3 = nonzero_constant_indices(large_eval_3);  
print(nonzero_large_eval_gen_3);

eval_4 = evaluate_variables(gens(minimal_ideal_4), [11, 12, 13]);
nonzero_eval_gen_4 = nonzero_constant_indices(eval_4);  # Shows that there is no non-zero constant generator when all of b[1], b[2], and b[3] are set to 0
print(nonzero_eval_gen_4);
large_eval_4 = evaluate_variables(gens(large_minimal_ideal_4), [11, 12, 13]);
nonzero_large_eval_gen_4 = nonzero_constant_indices(large_eval_4);  
print(nonzero_large_eval_gen_4);

# From the above evaluations, we can conclude that no union of charts covers one of the remaining ones

# We compute the radicals of the ideals inducing the affine charts
rad_chart_1 = radical(b_minimal_ideal_1);
rad_chart_2 = radical(b_minimal_ideal_2);
rad_chart_3 = radical(b_minimal_ideal_3);
rad_chart_4 = radical(b_minimal_ideal_4);

# We define the varieties corresponding to the charts to be able to compute their singular loci
X_1 = spec(C, rad_chart_1);
X_2 = spec(C, rad_chart_2);
X_3 = spec(C, rad_chart_3);
X_4 = spec(C, rad_chart_4);

sing_1 = singular_locus(X_1);
sing_2 = singular_locus(X_2);
sing_3 = singular_locus(X_3);
sing_4 = singular_locus(X_4);
sing_1_ideal = modulus(OO(sing_1[1]));
sing_2_ideal = modulus(OO(sing_2[1]));
sing_3_ideal = modulus(OO(sing_3[1]));
sing_4_ideal = modulus(OO(sing_4[1]));
