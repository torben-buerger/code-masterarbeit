using Oscar
Q_t, t = polynomial_ring(QQ, :t);
K, zeta = number_field(t^4 - t^2 + 1, "zeta");
i = zeta^3;      # imaginary unit
omega = zeta^4;  # primitive third root of unity

r = matrix(K, 2, 2, [1 0; 0 -1]);
r_1 = matrix(K, 2, 2, [(omega//2)*(-1-i) (omega//2)*(1-i); (omega//2)*(-1-i) (omega//2)*(-1+i)]);
G_6 = matrix_group([r, r_1]);

print(describe(G_6));

r_symp = block_diagonal_matrix([r, transpose(inv(r))]);
r_1_symp = block_diagonal_matrix([r_1, transpose(inv(r_1))]);
G_6_symp = matrix_group([r_symp, r_1_symp]);

# Compute the invariants of the symplectic reflection group
IR = invariant_ring(G_6_symp);
poly_invar = polynomial_ring(IR);
u = gens(poly_invar);
invars = fundamental_invariants(IR);
invars_ideal_1 = ideal(poly_invar, invars);

# Define the invariants as polynomials in the original variables
R, x = polynomial_ring(K, :x => (1:4));

invar_1 = x[1]*x[3] + x[2]*x[4];
invar_2 = x[1]^4 + (4*zeta^2 - 2)*x[1]^2*x[2]^2 + x[2]^4;
invar_3 = x[3]^4 + (-4*zeta^2 + 2)*x[3]^2*x[4]^2 + x[4]^4;
invar_4 = x[1]*x[2]^2*x[3]^3 - x[2]^3*x[3]^2*x[4] - x[1]^3*x[3]*x[4]^2 + x[1]^2*x[2]*x[4]^3;
invar_5 = x[1]^4*x[2]^2*x[3]^2 + (4//3*zeta^2 - 2//3)*x[1]^2*x[2]^4*x[3]^2 - 1//3*x[2]^6*x[3]^2 + (-4//3*zeta^2 + 2//3)*x[1]^5*x[2]*x[3]*x[4] - 4//3*x[1]^3*x[2]^3*x[3]*x[4] + (-4//3*zeta^2 + 2//3)*x[1]*x[2]^5*x[3]*x[4] - 1//3*x[1]^6*x[4]^2 + (4//3*zeta^2 - 2//3)*x[1]^4*x[2]^2*x[4]^2 + x[1]^2*x[2]^4*x[4]^2;
invar_6 = x[2]^2*x[3]^6 + (-4*zeta^2 + 2)*x[1]*x[2]*x[3]^5*x[4] - 3*x[1]^2*x[3]^4*x[4]^2 + (4*zeta^2 - 2)*x[2]^2*x[3]^4*x[4]^2 + 4*x[1]*x[2]*x[3]^3*x[4]^3 + (4*zeta^2 - 2)*x[1]^2*x[3]^2*x[4]^4 - 3*x[2]^2*x[3]^2*x[4]^4 + (-4*zeta^2 + 2)*x[1]*x[2]*x[3]*x[4]^5 + x[1]^2*x[4]^6;
invar_7 = x[1]^7*x[2]^2*x[3] + (2//3*zeta^2 - 1//3)*x[1]^5*x[2]^4*x[3] - x[1]^3*x[2]^6*x[3] + (-2//3*zeta^2 + 1//3)*x[1]*x[2]^8*x[3] + (-2//3*zeta^2 + 1//3)*x[1]^8*x[2]*x[4] - x[1]^6*x[2]^3*x[4] + (2//3*zeta^2 - 1//3)*x[1]^4*x[2]^5*x[4] + x[1]^2*x[2]^7*x[4];
invar_8 = x[2]*x[3]^8*x[4] + (-2*zeta^2 + 1)*x[1]*x[3]^7*x[4]^2 + (2*zeta^2 - 1)*x[2]*x[3]^6*x[4]^3 - x[1]*x[3]^5*x[4]^4 - x[2]*x[3]^4*x[4]^5 + (2*zeta^2 - 1)*x[1]*x[3]^3*x[4]^6 + (-2*zeta^2 + 1)*x[2]*x[3]^2*x[4]^7 + x[1]*x[3]*x[4]^8;
invar_9 = x[3]^12 - 33*x[3]^8*x[4]^4 - 33*x[3]^4*x[4]^8 + x[4]^12;
invar_10 = x[1]^12 - 33*x[1]^8*x[2]^4 - 33*x[1]^4*x[2]^8 + x[2]^12;

invars_ideal_2 = ideal(R, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, invar_9, invar_10]);

# Compute the relations between the invariants to obtain the presentation of the invariant ring
S, y = polynomial_ring(K, :y => (1:10));
pi = hom(S, R, [gens(invars_ideal_2)[i] for i in 1:10]);
relations_ideal = kernel(pi);
relations_basis = standard_basis(relations_ideal, ordering=negdegrevlex(S));
print(length(relations_basis));
relations_basis_ideal = ideal(S, relations_basis);

# Compute the conjugacy classes of reflections in G_6
classes = conjugacy_classes(G_6);
eigenvals = [eigenvalues(matrix(representative(c))) for c in classes];

ref_1 = representative(classes[2]);
ref_1_list = collect(classes[2]);
ref_2 = representative(classes[8]);
ref_2_list = collect(classes[8]);
ref_3 = representative(classes[13]);
ref_3_list = collect(classes[13]);

# Find the reflecting hyperplanes (= lines as we are in dimension 2) of the reflections by computing the eigenvectors. We have to take the transposed matrices because Oscar computes the action from the right
ref_1_eigvecs = [eigenspaces(K, transpose(matrix(m))) for m in ref_1_list];
ref_2_eigvecs = [eigenspaces(K, transpose(matrix(m))) for m in ref_2_list];
ref_3_eigvecs = [eigenspaces(K, transpose(matrix(m))) for m in ref_3_list];

# Extract the polynomials and vectors defining the hyperplanes from the eigenvectors given as column vectors corresponding to the eigenvalue 1
ref_poly_1 = ref_1_eigvecs[1][K(1)][2]*x[1] - ref_1_eigvecs[1][K(1)][1]*x[2];
ref_poly_2 = ref_2_eigvecs[1][K(1)][2]*x[1] - ref_2_eigvecs[1][K(1)][1]*x[2];
ref_poly_3 = ref_3_eigvecs[2][K(1)][2]*x[1] - ref_3_eigvecs[2][K(1)][1]*x[2];

ref_1_vec = matrix(K, 2, 1, [ref_1_eigvecs[1][K(1)][1]; ref_1_eigvecs[1][K(1)][2]]);
ref_2_vec = matrix(K, 2, 1, [ref_2_eigvecs[1][K(1)][1]; ref_2_eigvecs[1][K(1)][2]]);
ref_3_vec = matrix(K, 2, 1, [ref_3_eigvecs[2][K(1)][1]; ref_3_eigvecs[2][K(1)][2]]);
ref_vectors = [ref_1_vec, ref_2_vec, ref_3_vec];

# Compute the orbits of the action of G_6 on the reflecting hyperplanes by applying the group elements to the polynomials defining the hyperplanes
G_6_tr = matrix_group([transpose(matrix(m)) for m in gens(G_6)]);
G_6_transpose = collect(G_6_tr);
ref_1_orbit = [m*ref_1_vec for m in G_6_transpose];
ref_2_orbit = [m*ref_2_vec for m in G_6_transpose];
ref_3_orbit = [m*ref_3_vec for m in G_6_transpose];

# Compute unique representatives of all lines that form the orbits
function unique_lines(vectors)
    result_lines = [];
    push!(result_lines, vectors[1]);
    for v in vectors[2:end]
        is_multiple = false

        for w in result_lines
            # Check if v and w are scalar multiples
            if w[1,1] != 0
                lambda = v[1,1] / w[1,1]
                if v == lambda * w
                    is_multiple = true
                    break
                end
            elseif w[2,1] != 0
                lambda = v[2,1] / w[2,1]
                if v == lambda * w
                    is_multiple = true
                    break
                end
            end
        end

        if !is_multiple  # If v is not a scalar multiple of one of the previous vectors, add it to the list
            push!(result_lines, v)
        end
    end

    for v in result_lines
        if v[1,1] != 0
            v ./= v[1,1]  # Normalize the vector so that the first entry is 1
        elseif v[2,1] != 0
            v ./= v[2,1]  # Normalize the vector so that the second entry is 1
        end
    end
    return result_lines
end

ref_1_unique_lines = unique_lines(ref_1_orbit);
ref_2_unique_lines = unique_lines(ref_2_orbit);
ref_3_unique_lines = unique_lines(ref_3_orbit);

# Determine the inertia groups of the lines by checking which group elements fix the lines, the action in Oscar is from the right, so we have to transpose the matrix group
V = vector_space(K, 2);

# First, consider the orbit omega_1 corresponding to ref_1
v_11 = V(collect(ref_1_unique_lines[1][:, 1]));
v_12 = V(collect(ref_1_unique_lines[2][:, 1]));
v_13 = V(collect(ref_1_unique_lines[3][:, 1]));
v_14 = V(collect(ref_1_unique_lines[4][:, 1]));
v_15 = V(collect(ref_1_unique_lines[5][:, 1]));
v_16 = V(collect(ref_1_unique_lines[6][:, 1]));
stabilizer_11 = stabilizer(G_6_tr, v_11);
stabilizer_12 = stabilizer(G_6_tr, v_12);
stabilizer_13 = stabilizer(G_6_tr, v_13);
stabilizer_14 = stabilizer(G_6_tr, v_14);
stabilizer_15 = stabilizer(G_6_tr, v_15);
stabilizer_16 = stabilizer(G_6_tr, v_16);  # All of order 2, containing the corresponding reflection (which is its own inverse) and the identity

# Define the linear forms whose kernel are the reflecting hyperplanes and compute the corresponding polynomial delta_1 for the first orbit
ref_poly_11 = v_11[2]*x[1] - v_11[1]*x[2];
ref_poly_12 = v_12[2]*x[1] - v_12[1]*x[2];
ref_poly_13 = v_13[2]*x[1] - v_13[1]*x[2];
ref_poly_14 = v_14[2]*x[1] - v_14[1]*x[2];
ref_poly_15 = v_15[2]*x[1] - v_15[1]*x[2];
ref_poly_16 = v_16[2]*x[1] - v_16[1]*x[2];
delta_11 = ref_poly_11*ref_poly_12*ref_poly_13*ref_poly_14*ref_poly_15*ref_poly_16;  # Using the notation of the Bonnafé construction, this is the only polynomial for the first orbit since e_omega_1 = 2

# Next, consider the orbit omega_2 corresponding to ref_2 and ref_3
v_21 = V(collect(ref_2_unique_lines[1][:, 1]));
v_22 = V(collect(ref_2_unique_lines[2][:, 1]));
v_23 = V(collect(ref_2_unique_lines[3][:, 1]));
v_24 = V(collect(ref_2_unique_lines[4][:, 1]));
stabilizer_21 = stabilizer(G_6_tr, v_21);
stabilizer_22 = stabilizer(G_6_tr, v_22);
stabilizer_23 = stabilizer(G_6_tr, v_23);
stabilizer_24 = stabilizer(G_6_tr, v_24);  # All of order 3, containing the reflection, its inverse and the identity

# Define the linear forms whose kernel are the reflecting hyperplanes and compute the corresponding polynomials delta_1 and delta_2 for the first orbit
ref_poly_21 = v_21[2]*x[1] - v_21[1]*x[2];
ref_poly_22 = v_22[2]*x[1] - v_22[1]*x[2];
ref_poly_23 = v_23[2]*x[1] - v_23[1]*x[2];
ref_poly_24 = v_24[2]*x[1] - v_24[1]*x[2];
delta_21 = ref_poly_21*ref_poly_22*ref_poly_23*ref_poly_24;  # Using the notation of the Bonnafé construction, this is the first polynomial for the second orbit
delta_22 = ref_poly_21^2*ref_poly_22^2*ref_poly_23^2*ref_poly_24^2;  # This is the final polynomial for omega_2 since e_omega_2 = 3

# Save the results of the relations computation to save the time of recomputing them
relations_ideal_1 = ideal(S, [-3*y[1]^6 + (-20*zeta^2 + 10)*y[1]^3*y[4] + 3*y[1]^2*y[2]*y[3] + (2*zeta^2 - 1)*y[2]*y[6] + (6*zeta^2 - 3)*y[3]*y[5] + 9*y[4]^2,
  y[1]^5*y[2] + (-6*zeta^2 + 3)*y[1]^3*y[5] + (6*zeta^2 - 3)*y[1]^2*y[2]*y[4] - y[1]*y[2]^2*y[3] + (4*zeta^2 - 2)*y[3]*y[7] + 3*y[4]*y[5],
  y[1]^4*y[2]^2 + (4*zeta^2 - 2)*y[1]^3*y[7] + (-6*zeta^2 + 3)*y[1]^2*y[2]*y[5] + (6*zeta^2 - 3)*y[1]*y[2]^2*y[4] - y[2]^3*y[3] - 18*y[4]*y[7] + 9*y[5]^2,
  y[1]^5*y[3] + (-2*zeta^2 + 1)*y[1]^3*y[6] + (6*zeta^2 - 3)*y[1]^2*y[3]*y[4] - y[1]*y[2]*y[3]^2 + 2*y[2]*y[8] + y[4]*y[6],
  (16*zeta^2 - 8)*y[1]^5*y[4] - 3*y[1]^4*y[2]*y[3] + (2*zeta^2 - 1)*y[1]^2*y[2]*y[6] + (6*zeta^2 - 3)*y[1]^2*y[3]*y[5] - 72*y[1]^2*y[4]^2 + (-36*zeta^2 + 18)*y[1]*y[2]*y[3]*y[4] + 3*y[2]^2*y[3]^2 + 9*y[5]*y[6],
  y[1]^4*y[3]^2 + 2*y[1]^3*y[8] + (-2*zeta^2 + 1)*y[1]^2*y[3]*y[6] + (6*zeta^2 - 3)*y[1]*y[3]^2*y[4] - y[2]*y[3]^3 + (12*zeta^2 - 6)*y[4]*y[8] + y[6]^2,
  (-36*zeta^2 + 18)*y[1]^3*y[7] + (18*zeta^2 - 9)*y[1]^2*y[2]*y[5] + (-6*zeta^2 + 3)*y[1]*y[2]^2*y[4] + 2*y[2]^3*y[3] - 2*y[3]*y[10] + 18*y[4]*y[7],
  2*y[1]^3*y[2]^3 - 2*y[1]^3*y[10] + (-36*zeta^2 + 18)*y[1]^2*y[2]*y[7] + (18*zeta^2 - 9)*y[1]*y[2]^2*y[5] + (6*zeta^2 - 3)*y[2]^3*y[4] + (-12*zeta^2 + 6)*y[4]*y[10] + 54*y[5]*y[7],
  (24*zeta^2 - 12)*y[1]^5*y[5] + (-12*zeta^2 + 6)*y[1]^2*y[3]*y[7] - 36*y[1]^2*y[4]*y[5] + (-2*zeta^2 + 1)*y[1]*y[2]^2*y[6] + (-12*zeta^2 + 6)*y[1]*y[2]*y[3]*y[5] + (2*zeta^2 - 1)*y[2]^2*y[3]*y[4] + 6*y[6]*y[7],
  (2*zeta^2 - 1)*y[2]^3*y[5] + (-2*zeta^2 + 1)*y[5]*y[10] + 18*y[7]^2,
  (-36*zeta^2 + 18)*y[1]^3*y[8] - 9*y[1]^2*y[3]*y[6] + 9*y[1]*y[3]^2*y[4] + (4*zeta^2 - 2)*y[2]*y[3]^3 + (-4*zeta^2 + 2)*y[2]*y[9] + 18*y[4]*y[8],
  -4*y[1]^5*y[6] + (-4*zeta^2 + 2)*y[1]^2*y[2]*y[8] + (-8*zeta^2 + 4)*y[1]^2*y[4]*y[6] + 2*y[1]*y[2]*y[3]*y[6] + 3*y[1]*y[3]^2*y[5] - y[2]*y[3]^2*y[4] + 6*y[5]*y[8],
  (4*zeta^2 - 2)*y[1]^3*y[3]^3 + (-4*zeta^2 + 2)*y[1]^3*y[9] + (-36*zeta^2 + 18)*y[1]^2*y[3]*y[8] - 9*y[1]*y[3]^2*y[6] - 9*y[3]^3*y[4] + 18*y[4]*y[9] + 18*y[6]*y[8],
  (-6*zeta^2 + 3)*y[1]^2*y[5]*y[6] + (2*zeta^2 - 1)*y[1]*y[2]*y[4]*y[6] + (6*zeta^2 - 3)*y[1]*y[3]*y[4]*y[5] + (-2*zeta^2 + 1)*y[2]*y[3]*y[4]^2 + 12*y[7]*y[8],
  (-2*zeta^2 + 1)*y[3]^3*y[6] + (2*zeta^2 - 1)*y[6]*y[9] + 18*y[8]^2,
  (2*zeta^2 - 1)*y[1]^2*y[6]^2 + (-4*zeta^2 + 2)*y[1]*y[3]*y[4]*y[6] - 2*y[3]^3*y[5] + (2*zeta^2 - 1)*y[3]^2*y[4]^2 + 2*y[5]*y[9],
  72*y[1]^2*y[5]*y[8] + (-6*zeta^2 + 3)*y[1]*y[2]*y[6]^2 + (-18*zeta^2 + 9)*y[1]*y[3]*y[5]*y[6] - 24*y[1]*y[4]^2*y[6] + (6*zeta^2 - 3)*y[2]*y[3]*y[4]*y[6] - 4*y[3]^3*y[7] + (18*zeta^2 - 9)*y[3]^2*y[4]*y[5]+ 24*y[3]*y[4]^3 + 4*y[7]*y[9],
  (54*zeta^2 - 27)*y[1]^2*y[5]^2 + (-36*zeta^2 + 18)*y[1]*y[2]*y[4]*y[5] - 2*y[2]^3*y[6] + (6*zeta^2 - 3)*y[2]^2*y[4]^2 + 2*y[6]*y[10],
  -72*y[1]^2*y[6]*y[7] + 27*y[1]*y[2]*y[5]*y[6] + 81*y[1]*y[3]*y[5]^2 + (-144*zeta^2 + 72)*y[1]*y[4]^2*y[5] - 4*y[2]^3*y[8] - 9*y[2]^2*y[4]*y[6] - 27*y[2]*y[3]*y[4]*y[5] + (48*zeta^2 - 24)*y[2]*y[4]^3 + 4*y[8]*y[10],
  -3290101285742342628191948215568900615539620384719577*y[1]^4*y[5]*y[6] + 30232015654944220871070857400173055270528000*y[1]^3*y[2]^2*y[8] - 367565213137716042200644235006599469144736404400429*y[1]^3*y[2]*y[4]*y[6] 
  + (12092806261977688348428342960069222108211200*zeta^2 - 6046403130988844174214171480034611054105600)*y[1]^3*y[3]^2*y[7] + (-12092806261977688348428342960069222108211200*zeta^2 + 6046403130988844174214171480034611054105600)*y[1]^2*y[2]^2*y[3]*y[6] 
  + (-870682050862393561086840693124983991791206400*zeta^2 + 435341025431196780543420346562491995895603200)*y[1]^2*y[7]*y[8] + (-3923530756510833051631564619085232598474160487405968*zeta^2 
  + 1961765378255416525815782309542616299237080243702984)*y[1]*y[2]*y[5]*y[8] + 228637438148936169846427640834555114170169574692820*y[1]*y[3]*y[6]*y[7] + (-2193400421820536320931184933625587181201084360876518*zeta^2 + 1096700210910268160465592466812793590600542180438259)*y[1]*y[4]*y[5]*y[6] 
  + (-438612412075062045731395513436405042383092866794720*zeta^2 + 219306206037531022865697756718202521191546433397360)*y[2]^2*y[4]*y[8] + 38605887453810897682776120840324117262236355513989*y[2]^2*y[6]^2 + 1098199285706500867598067798074718076054956886430616*y[2]*y[3]*y[5]*y[6] 
  + (-450941727268498166099790143720720325889823824913294*zeta^2 + 225470863634249083049895071860360162944911912456647)*y[2]*y[4]^2*y[6] - 2015467710329614724738057160011537018035200*y[3]^3*y[10] - 654921048974828229261183332191220675988471082344966*y[3]^2*y[4]*y[7] 
  + 2947144870035204523649218306661237172804743459665947*y[3]^2*y[5]^2 + (-5239369262480676696483027744370458532891760449966128*zeta^2 + 2619684631240338348241513872185229266445880224983064)*y[3]*y[4]^2*y[5] 
  - 870682050862393561086840693124983991791206400*y[4]^4 + 2015467710329614724738057160011537018035200*y[9]*y[10]]);

# Define the quotient ring of S by the relations ideal which is isomorphic to the invariant ring of G_5_symp
S_quo, quo_map = quo(S, relations_ideal_1);

# Define the ideals corresponding to the polynomials that arise from the reflecting hyperplanes in the Bonnafé construction
delta_11_ideal = ideal(R, delta_11);
delta_21_ideal = ideal(R, delta_21);
delta_22_ideal = ideal(R, delta_22);
delta_11_21_ideal = ideal(R, delta_11*delta_21);

# Investigate the first polynomial delta_11 and compute the variety along which the blow-up shall be performed
W_delta_11 = preimage(pi, delta_11_ideal);
basis_W_delta_11 = standard_basis(W_delta_11, ordering=negdegrevlex(S));
quo_W_delta_11 = quo_map(W_delta_11);
quo_basis_W_delta_11 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_11]);
print(length(basis_W_delta_11));

nonzero_quo_W_delta_11 = gens(quo_basis_W_delta_11)[[1:7; 9; 12; 16:18]];
test_delta_11 = ideal(S_quo, nonzero_quo_W_delta_11[1:4]);
print(test_delta_11 == quo_W_delta_11);  # Returns true, therefore, the ideal image of W_delta_11 in the quotient ring is generated by the first 4 nonzero elements of the basis
rad_test_delta_11 = radical(test_delta_11);
print(rad_test_delta_11 == test_delta_11);  # Retuns true and thus shows that the ideal is radical, so blowing up along this ideal is equivalent to blowing up the corresponding subvariety

# List the generators of the ideal W_delta_11_prime which satisfies W_delta_11 = W_delta_11_prime + relations_ideal
w_11_1 = y[7];
w_11_2 = -y[2]^3 + y[10];
w_11_3 = -3*y[1]*y[5] + y[2]*y[4];
w_11_4 = 4*y[1]^3*y[4] - y[2]*y[6] - 3*y[3]*y[5] + (8*zeta^2-4)*y[4]^2;

W_delta_11_prime = ideal(S, [w_11_1, w_11_2, w_11_3, w_11_4]);
print(W_delta_11 == W_delta_11_prime + relations_ideal);  # Returns true


# Investigate the third polynomial delta_21 and compute the variety along which the blow-up shall be performed
W_delta_21 = preimage(pi, delta_21_ideal);
basis_W_delta_21 = standard_basis(W_delta_21, ordering=negdegrevlex(S));
quo_W_delta_21 = quo_map(W_delta_21);
quo_basis_W_delta_21 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_21]);
print(length(basis_W_delta_21));

nonzero_quo_W_delta_21 = gens(quo_basis_W_delta_21)[1:8];
test_delta_21 = ideal(S_quo, nonzero_quo_W_delta_21[1:4]);
print(test_delta_21 == quo_W_delta_21);  # Returns true, therefore, the ideal image of W_delta_21 in the quotient ring is generated by the first 4 nonzero elements of the basis
rad_test_delta_21 = radical(test_delta_21);
print(rad_test_delta_21 == test_delta_21);  # Retuns true and thus shows that the ideal is radical, so blowing up along this ideal is equivalent to blowing up the corresponding subvariety

# List the generators of the ideal W_delta_21_prime which satisfies W_delta_21 = W_delta_21_prime + relations_ideal
w_21_1 = y[1]^2*y[2] + (-6*zeta^2 + 3)*y[5];
w_21_2 = y[1]*y[2]^2 + (-12*zeta^2 + 6)*y[7];
w_21_3 = -y[2]^3 + 2*y[10];
w_21_4 = (-8*zeta^2 + 4)*y[1]*y[4] + y[2]*y[3];

W_delta_21_prime = ideal(S, [w_21_1, w_21_2, w_21_3, w_21_4]);
print(W_delta_21 == W_delta_21_prime + relations_ideal);  # Returns true


# Investigate the final polynomial delta_22 and compute the variety along which the blow-up shall be performed
W_delta_22 = preimage(pi, delta_22_ideal);
basis_W_delta_22 = standard_basis(W_delta_22, ordering=negdegrevlex(S));
quo_W_delta_22 = quo_map(W_delta_22);
quo_basis_W_delta_22 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_22]);
print(length(basis_W_delta_22));

nonzero_quo_W_delta_22 = gens(quo_basis_W_delta_22)[[1; 3:5; 9:13; 17; 20:21;]];
test_delta_22 = ideal(S_quo, nonzero_quo_W_delta_22[[1:3; 12]]);
print(test_delta_22 == quo_W_delta_22);  # Returns true, therefore, the ideal image of W_delta_22 in the quotient ring is generated by 4 elements of the basis
rad_test_delta_22 = radical(test_delta_22);
print(rad_test_delta_22 == test_delta_22);  

# List the generators of the ideal W_delta_21_prime which satisfies W_delta_21 = W_delta_21_prime + relations_ideal
w_22_1 = -y[2]^3 + 2*y[10];
w_22_2 = (-2*zeta^2 + 1)*y[1]^2*y[2]^2 - 36*y[1]*y[7] + 9*y[2]*y[5];
w_22_3 = (4*zeta^2 - 2)*y[1]^2*y[2]*y[4] - y[1]*y[2]^2*y[3] + (12*zeta^2 - 6)*y[3]*y[7] + 18*y[4]*y[5];
w_22_4 = -48*y[1]^2*y[4]^2 + (-16*zeta^2 + 8)*y[1]*y[2]*y[3]*y[4] + y[2]^2*y[3]^2;

W_delta_22_prime = ideal(S, [w_22_1, w_22_2, w_22_3, w_22_4]);
print(W_delta_22 == W_delta_22_prime + relations_ideal);  # Returns true
