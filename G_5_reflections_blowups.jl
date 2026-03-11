using Oscar
Q_t, t = polynomial_ring(QQ, :t);
K, zeta = number_field(t^4 - t^2 + 1, "zeta");
i = zeta^3;      # imaginary unit
sigma = zeta^4;  # primitive third root of unity

# Matrix representations of generators of the group G_5 as subgroup of GL_2(CC) using the generators given by Leher--Taylor
r = matrix(K, 2, 2, [1 0; 0 -1]);
r_1 = (sigma//2)*matrix([-1-i 1-i; -1-i -1+i]);
r_2 = (sigma//2)*matrix([-1+i -1+i; 1+i -1-i]);
r_2_conj = r*r_2*r;

G_5 = matrix_group([r_1, r_2_conj]);
describe(G_5);

# Find the conjugacy classes of G_5 consisting of reflections
classes = conjugacy_classes(G_5);
eigenvals = [eigenvalues(matrix(representative(c))) for c in classes];

ref_1 = representative(classes[11]);
ref_1_list = collect(classes[11]);
ref_2 = representative(classes[12]);
ref_2_list = collect(classes[12]);
ref_3 = representative(classes[16]);
ref_3_list = collect(classes[16]);
ref_4 = representative(classes[17]);
ref_4_list = collect(classes[17]);

# Find the reflecting hyperplanes (= lines as we are in dimension 2) of the reflections by computing the eigenvectors 
ref_1_eigvecs = [eigenspaces(K, matrix(m)) for m in ref_1_list];
ref_2_eigvecs = [eigenspaces(K, matrix(m)) for m in ref_2_list];
ref_3_eigvecs = [eigenspaces(K, matrix(m)) for m in ref_3_list];
ref_4_eigvecs = [eigenspaces(K, matrix(m)) for m in ref_4_list];

# Extract the polynomials and vectors defining the hyperplanes from the eigenvectors corresponding to the eigenvalue 1
R, x = polynomial_ring(K, :x => (1:2));
ref_poly_1 = ref_1_eigvecs[1][K(1)][2]*x[1] - ref_1_eigvecs[1][K(1)][1]*x[2];
ref_poly_2 = ref_2_eigvecs[1][K(1)][2]*x[1] - ref_2_eigvecs[1][K(1)][1]*x[2];
ref_poly_3 = ref_3_eigvecs[2][K(1)][2]*x[1] - ref_3_eigvecs[2][K(1)][1]*x[2];
ref_poly_4 = ref_4_eigvecs[1][K(1)][2]*x[1] - ref_4_eigvecs[1][K(1)][1]*x[2];

ref_1_vec = matrix(K, 2, 1, [ref_1_eigvecs[1][K(1)][1]; ref_1_eigvecs[1][K(1)][2]]);
ref_2_vec = matrix(K, 2, 1, [ref_2_eigvecs[1][K(1)][1]; ref_2_eigvecs[1][K(1)][2]]);
ref_3_vec = matrix(K, 2, 1, [ref_3_eigvecs[2][K(1)][1]; ref_3_eigvecs[2][K(1)][2]]);
ref_4_vec = matrix(K, 2, 1, [ref_4_eigvecs[1][K(1)][1]; ref_4_eigvecs[1][K(1)][2]]);
ref_vectors = [ref_1_vec, ref_2_vec, ref_3_vec, ref_4_vec];

# Compute the orbits of the action of G_5 on the reflecting hyperplanes by applying the group elements to the polynomials defining the hyperplanes
G_5_list = collect(G_5);
ref_1_orbit = [matrix(m)*ref_1_vec for m in G_5_list];
ref_2_orbit = [matrix(m)*ref_2_vec for m in G_5_list];
ref_3_orbit = [matrix(m)*ref_3_vec for m in G_5_list];
ref_4_orbit = [matrix(m)*ref_4_vec for m in G_5_list];

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
ref_4_unique_lines = unique_lines(ref_4_orbit);  # The orbits of ref_1 and ref_4 are the same, as well as those of ref_2 and ref_3

# Determine the inertia groups of the lines by checking which group elements fix the lines, the action in Oscar is from the right, so we have to transpose the matrix group
G_5_transpose = matrix_group([transpose(matrix(m)) for m in gens(G_5)]);
V = vector_space(K, 2);
# First, consider the orbit omega_1 corresponding to ref_1 and ref_4
v_11 = V(collect(ref_1_unique_lines[1][:, 1]));
v_12 = V(collect(ref_1_unique_lines[2][:, 1]));
v_13 = V(collect(ref_1_unique_lines[3][:, 1]));
v_14 = V(collect(ref_1_unique_lines[4][:, 1]));
stabilizer_11 = stabilizer(G_5_transpose, v_11);
stabilizer_12 = stabilizer(G_5_transpose, v_12);
stabilizer_13 = stabilizer(G_5_transpose, v_13);
stabilizer_14 = stabilizer(G_5_transpose, v_14);  # All of order 3, containing the reflection, its inverse and the identity
ref_poly_11 = v_11[2]*x[1] - v_11[1]*x[2];
ref_poly_12 = v_12[2]*x[1] - v_12[1]*x[2];
ref_poly_13 = v_13[2]*x[1] - v_13[1]*x[2];
ref_poly_14 = v_14[2]*x[1] - v_14[1]*x[2];
delta_11 = ref_poly_11*ref_poly_12*ref_poly_13*ref_poly_14;  # Using the notation of the Bonnafé construction, this is the polynomial for the first orbit
delta_12 = ref_poly_11^2*ref_poly_12^2*ref_poly_13^2*ref_poly_14^2;  # This is the final polynomial for omega_1 since e_omega_1 = 3

# Next, consider th orbit omega_2 corresponding to ref_2 and ref_3
v_21 = V(collect(ref_2_unique_lines[1][:, 1]));
v_22 = V(collect(ref_2_unique_lines[2][:, 1]));
v_23 = V(collect(ref_2_unique_lines[3][:, 1]));
v_24 = V(collect(ref_2_unique_lines[4][:, 1]));
stabilizer_21 = stabilizer(G_5_transpose, v_21);
stabilizer_22 = stabilizer(G_5_transpose, v_22);
stabilizer_23 = stabilizer(G_5_transpose, v_23);
stabilizer_24 = stabilizer(G_5_transpose, v_24);  # All of order 3, containing the reflection, its inverse and the identity
ref_poly_21 = v_21[2]*x[1] - v_21[1]*x[2];
ref_poly_22 = v_22[2]*x[1] - v_22[1]*x[2];
ref_poly_23 = v_23[2]*x[1] - v_23[1]*x[2];
ref_poly_24 = v_24[2]*x[1] - v_24[1]*x[2];
delta_21 = ref_poly_21*ref_poly_22*ref_poly_23*ref_poly_24;  # Using the notation of the Bonnafé construction, this is the polynomial for the second orbit
delta_22 = ref_poly_21^2*ref_poly_22^2*ref_poly_23^2*ref_poly_24^2;  # This is the final polynomial for omega_2 since e_omega_2 = 3

# Compute the coordinate ring of the variety (V+V^*)/G_5 by finding the invariant polynomials
r_1_symp = block_diagonal_matrix([r_1, transpose(inv(r_1))]);
r_2_conj_symp = block_diagonal_matrix([r_2_conj, transpose(inv(r_2_conj))]);
G_5_symp = matrix_group([r_1_symp, r_2_conj_symp]);
IR = invariant_ring(G_5_symp);
poly_invar = polynomial_ring(IR);
u = gens(poly_invar);
invars = fundamental_invariants(IR);
invars_ideal_1 = ideal(poly_invar, invars);

# Compute the relations of the invariants and therefore a presentation of the invariant ring as quotient of a polynomial ring
B, y = polynomial_ring(K, :y => (1:12));
f = hom(B, poly_invar, invars);
# The kernel of f is the ideal of relations that has to be factored out to present the invariant ring
relations_ideal = kernel(f);  # Computation runs out of memory after three hours
