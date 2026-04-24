# Here, we consider G_5 using the presentation given by Berry to be able to using the relations given in the appendix of his paper
using Oscar
Q_t, t = polynomial_ring(QQ, :t);
K, zeta = number_field(t^8 - t^4 + 1, "zeta");
i = zeta^6;      # imaginary unit
sigma = zeta^8;  # primitive third root of unity
q2 = zeta^3 + zeta^(-3); # sqrt(2)
q3 = zeta^8 - zeta^(-8); # sqrt(3)

# Matrix representations of generators of the group G_5 as subgroup of GL_2(CC)
a = matrix(K, 2, 2, [i 0; 0 -i]);
b = matrix(K, 2, 2, [0 1; -1 0]);
c = matrix(K, 2, 2, [zeta//q2 (-zeta)//q2; (i*zeta)//q2 (i*zeta)//q2]);
d = matrix(K, 2, 2, [sigma, 0, 0, sigma]);

c_conj = a*c*d^2

G_5 = matrix_group([a, b, c, c_conj]);
describe(G_5);

# Compute the corresponding symplectic reflection group
a_symp = block_diagonal_matrix([a, transpose(inv(a))]);
b_symp = block_diagonal_matrix([b, transpose(inv(b))]);
c_symp = block_diagonal_matrix([c, transpose(inv(c))]);
d_symp = block_diagonal_matrix([d, transpose(inv(d))]);
c_conj_symp = block_diagonal_matrix([c_conj, transpose(inv(c_conj))]);
G_5_symp = matrix_group([a_symp, b_symp, c_symp, c_conj_symp]);

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

# Extract the polynomials and vectors defining the hyperplanes from the eigenvectors given as column vectors corresponding to the eigenvalue 1
R, x = polynomial_ring(K, :x => (1:4));
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

# Define the linear forms whose kernel are the reflecting hyperplanes and compute the corresponding polynomials delta_1 and delta_2 for the first orbit
ref_poly_11 = v_11[2]*x[1] - v_11[1]*x[2];
ref_poly_12 = v_12[2]*x[1] - v_12[1]*x[2];
ref_poly_13 = v_13[2]*x[1] - v_13[1]*x[2];
ref_poly_14 = v_14[2]*x[1] - v_14[1]*x[2];
delta_11 = ref_poly_11*ref_poly_12*ref_poly_13*ref_poly_14;  # Using the notation of the Bonnafé construction, this is the polynomial for the first orbit
delta_12 = ref_poly_11^2*ref_poly_12^2*ref_poly_13^2*ref_poly_14^2;  # This is the final polynomial for omega_1 since e_omega_1 = 3

# Next, consider the orbit omega_2 corresponding to ref_2 and ref_3
v_21 = V(collect(ref_2_unique_lines[1][:, 1]));
v_22 = V(collect(ref_2_unique_lines[2][:, 1]));
v_23 = V(collect(ref_2_unique_lines[3][:, 1]));
v_24 = V(collect(ref_2_unique_lines[4][:, 1]));
stabilizer_21 = stabilizer(G_5_transpose, v_21);
stabilizer_22 = stabilizer(G_5_transpose, v_22);
stabilizer_23 = stabilizer(G_5_transpose, v_23);
stabilizer_24 = stabilizer(G_5_transpose, v_24);  # All of order 3, containing the reflection, its inverse and the identity

# Define the linear forms whose kernel are the reflecting hyperplanes and compute the corresponding polynomials delta_1 and delta_2 for the first orbit
ref_poly_21 = v_21[2]*x[1] - v_21[1]*x[2];
ref_poly_22 = v_22[2]*x[1] - v_22[1]*x[2];
ref_poly_23 = v_23[2]*x[1] - v_23[1]*x[2];
ref_poly_24 = v_24[2]*x[1] - v_24[1]*x[2];
delta_21 = ref_poly_21*ref_poly_22*ref_poly_23*ref_poly_24;  # Using the notation of the Bonnafé construction, this is the polynomial for the second orbit
delta_22 = ref_poly_21^2*ref_poly_22^2*ref_poly_23^2*ref_poly_24^2;  # This is the final polynomial for omega_2 since e_omega_2 = 3

# Define the coordinate ring of the variety (V+V^*)/G_5_symp by taking the relations of the invariant ring of G_5_symp computed by Berry, see G_5_construction.jl
# First, we define the invariants of the symplectic action of G_5, see G_5_construction.jl for details
# Eigenvectors of c and c_conj for the eigenvalue 1 according to Berry
F11 = x[1]*x[3] + x[2]*x[4];
F40 = x[1]^4 + x[2]^4 - 2*q3*x[1]^2*x[2]^2;
F31 = x[1]^3*x[4] - q3*x[1]^2*x[2]*x[3] + q3*x[1]*x[2]^2*x[4] - x[2]^3*x[3];
G22 = x[1]^2*x[3]^2 + q3*x[2]^2*x[3]^2 + q3*x[1]^2*x[4]^2 + x[2]^2*x[4]^2 - 4*x[1]*x[2]*x[3]*x[4];
F13 = x[1]*x[4]^3 + q3*x[2]*x[3]*x[4]^2 - q3*x[1]*x[3]^2*x[4] - x[2]*x[3]^3;
F04 = x[3]^4 + x[4]^4 + 2*q3*x[3]^2*x[4]^2;
F60 = x[1]^5*x[2] - x[1]*x[2]^5;
G51 = x[1]^5*x[3] - 5*x[1]^4*x[2]*x[4] - 5*x[1]*x[2]^4*x[3] + x[2]^5*x[4];
G42 = x[1]^4*x[3]*x[4] - 2*x[1]^3*x[2]*x[4]^2 + 2*x[1]*x[2]^3*x[3]^2 - x[2]^4*x[3]*x[4];
F33 = x[1]^3*x[3]*x[4]^2 - x[1]^2*x[2]*x[4]^3 - x[1]*x[2]^2*x[3]^3 + x[2]^3*x[3]^2*x[4];
G24 = 2*x[1]^2*x[3]*x[4]^3 + x[1]*x[2]*x[3]^4 - x[1]*x[2]*x[4]^4 - 2*x[2]^2*x[3]^3*x[4];
G15 = x[1]*x[3]^5 - 5*x[1]*x[3]*x[4]^4 - 5*x[2]*x[3]^4*x[4] + x[2]*x[4]^5;
F06 = x[3]^5*x[4] - x[3]*x[4]^5;

f_conj = hom(K, K, -zeta^7+zeta^3);
function conj_poly(poly)
    return map_coefficients(f_conj, poly)
end

F40_conj = conj_poly(F40);
F31_conj = conj_poly(F31);
G22_conj = conj_poly(G22);
F13_conj = conj_poly(F13);
F04_conj = conj_poly(F04);

# Invariants as given by Berry
B1 = -F40*F31 - F40_conj*F31_conj;
B0 = F31*F13 + F31_conj*F13_conj;
Bneg1 = -F13*F04 - F13_conj*F04_conj;
C2 = F40^3 + F40_conj^3;
C1 = -F31^3 - F31_conj^3;
C0 = F40*F13^2 + F40_conj*F13_conj^2 + 4*F11^3*F33;
Cneg1 = -F13^3 - F13_conj^3;
Cneg2 = F04^3 + F04_conj^3;
h = F11;
A1 = 6*F60;
A0 = 3*F33;
Aneg1 = 6*F06;

invars_ideal = ideal(R, [h, A1, A0, Aneg1, B1, B0, Bneg1, C2, C1, C0, Cneg1, Cneg2]);

# Finally, we define the polynomial ring of which the invariant ring is the quotient by dividing out the relations ideal, for simplicity, we use y[1], ..., y[12] as variables instead of (y_h, y_A1, y_A0, y_Aneg1, y_B1, y_B0, y_Bneg1, y_C2, y_C1, y_C0, y_Cneg1, y_Cneg2)
S, y = polynomial_ring(K, :y => (1:12));

# The map pi_1 maps the variables y[1], ..., y[12] to the invariants h, A1, A0, Aneg1, B1, B0, Bneg1, C2, C1, C0, Cneg1, Cneg2
pi_1 = hom(S, R, [h, A1, A0, Aneg1, B1, B0, Bneg1, C2, C1, C0, Cneg1, Cneg2]);

# We make use of the relations given by Berry
relations_ideal = ideal(S, [
    2*y[2]*y[4] - 2*y[3]^2 + 9*y[1]^2*y[6],
    y[2]*y[6] - y[3]*y[5] + 3*y[1]*y[9],
    y[4]*y[6] - y[3]*y[7] - 3*y[1]*y[11],
    y[2]*y[7] - y[3]*y[6] + 3*y[1]*y[10] + 2*y[1]^4*y[3],
    y[4]*y[5] - y[3]*y[6] - 3*y[1]*y[10] + 2*y[1]^4*y[3],
    2*y[5]*y[7] - 2*y[6]^2 - 6*y[1]^2*y[3]^2 - y[1]^4*y[6],
    6*y[2]*y[9] - 6*y[3]*y[8] + 9*y[1]*y[5]^2 - y[1]^3*y[2]^2,
    6*y[4]*y[11] - 6*y[3]*y[12] - 9*y[1]*y[7]^2 + y[1]^3*y[4]^2,
    6*y[2]*y[10] - 6*y[3]*y[9] + 9*y[1]*y[5]*y[6] - 4*y[1]^3*y[2]*y[3],
    6*y[4]*y[10] - 6*y[3]*y[11] - 9*y[1]*y[7]*y[6] + 4*y[1]^3*y[4]*y[3],
    6*y[2]*y[11] - 6*y[3]*y[10] + 9*y[1]*y[6]^2 - 4*y[1]^3*y[3]^2,
    6*y[4]*y[9] - 6*y[3]*y[10] - 9*y[1]*y[6]^2 + 4*y[1]^3*y[3]^2,
    6*y[2]*y[12] - 6*y[3]*y[11] + 9*y[1]*y[6]*y[7] - 28*y[1]^3*y[3]*y[4] - 36*y[1]^5*y[7],
    6*y[4]*y[8] - 6*y[3]*y[9] - 9*y[1]*y[6]*y[5] + 28*y[1]^3*y[3]*y[2] + 36*y[1]^5*y[5],
    3*y[5]*y[9] - 3*y[6]*y[8] + 3*y[1]*y[2]^2*y[3] + 4*y[1]^3*y[2]*y[5],
    3*y[7]*y[11] - 3*y[6]*y[12] - 3*y[1]*y[4]^2*y[3] - 4*y[1]^3*y[4]*y[7],
    6*y[5]*y[10] - 6*y[6]*y[9] + 6*y[1]*y[2]*y[3]^2 + 4*y[1]^3*y[3]*y[5] + y[1]^3*y[2]*y[6],
    6*y[7]*y[10] - 6*y[6]*y[11] - 6*y[1]*y[4]*y[3]^2 - 4*y[1]^3*y[3]*y[7] - y[1]^3*y[4]*y[6],
    6*y[5]*y[11] - 6*y[6]*y[10] + 6*y[1]*y[3]^3 + 5*y[1]^3*y[3]*y[6],
    6*y[7]*y[9] - 6*y[6]*y[10] - 6*y[1]*y[3]^3 - 5*y[1]^3*y[3]*y[6],
    6*y[5]*y[12] - 6*y[6]*y[11] + 6*y[1]*y[3]^2*y[4] + 35*y[1]^3*y[4]*y[6] - 84*y[1]^4*y[11] + 4*y[1]^7*y[4],
    6*y[7]*y[8] - 6*y[6]*y[9] - 6*y[1]*y[3]^2*y[2] - 35*y[1]^3*y[2]*y[6] - 84*y[1]^4*y[9] - 4*y[1]^7*y[2],
    36*y[8]*y[10] - 36*y[9]^2 + 54*y[1]^2*y[2]*y[3]*y[5] - 18*y[1]^3*y[8]*y[3] + 63*y[1]^4*y[5]^2 + y[1]^6*y[2]^2,
    36*y[12]*y[10] - 36*y[11]^2 + 54*y[1]^2*y[4]*y[3]*y[7] + 18*y[1]^3*y[12]*y[3] + 63*y[1]^4*y[7]^2 + y[1]^6*y[4]^2,
    9*y[8]*y[11] - 9*y[9]*y[10] + 27*y[1]^2*y[2]*y[3]*y[6] + 36*y[1]^3*y[9]*y[3] + 18*y[1]^4*y[5]*y[6] + 2*y[1]^6*y[2]*y[3],
    9*y[12]*y[9] - 9*y[11]*y[10] + 27*y[1]^2*y[4]*y[3]*y[6] - 36*y[1]^3*y[11]*y[3] + 18*y[1]^4*y[7]*y[6] + 2*y[1]^6*y[4]*y[3],
    2*y[8]*y[12] - 2*y[9]*y[11] + 9*y[1]^2*y[2]*y[4]*y[6] + 48*y[1]^4*y[6]^2 + 32*y[1]^6*y[2]*y[4] + 132*y[1]^8*y[6] - 8*y[1]^12,
    18*y[9]*y[11] - 18*y[10]^2 + 27*y[1]^2*y[2]*y[4]*y[6] + 126*y[1]^4*y[6]^2 + 8*y[1]^6*y[3]^2,
    6*y[8]*y[9] - 3*y[5]^3 + 2*y[2]^3*y[3] + 3*y[1]^2*y[2]^2*y[5],
    6*y[12]*y[11] - 3*y[7]^3 + 2*y[4]^3*y[3] + 3*y[1]^2*y[4]^2*y[7],
    6*y[9]^2 - 3*y[5]^2*y[6] + 2*y[2]^3*y[4] + 12*y[1]^2*y[2]*y[3]*y[5] - 28*y[1]^3*y[2]*y[9],
    6*y[11]^2 - 3*y[7]^2*y[6] + 2*y[4]^3*y[2] + 12*y[1]^2*y[4]*y[3]*y[7] + 28*y[1]^3*y[4]*y[11],
    6*y[9]*y[10] - 3*y[5]*y[6]^2 + 2*y[2]*y[3]^3 + 3*y[1]^2*y[2]*y[3]*y[6] + 4*y[1]^3*y[3]*y[9],
    6*y[11]*y[10] - 3*y[7]*y[6]^2 + 2*y[4]*y[3]^3 + 3*y[1]^2*y[4]*y[3]*y[6] - 4*y[1]^3*y[3]*y[11],
    18*y[10]^2 - 9*y[5]*y[6]*y[7] + 6*y[3]^4 + 9*y[1]^2*y[3]^2*y[6] - 8*y[1]^6*y[3]^2
]);

# Define the quotient ring of S by the relations ideal which is isomorphic to the invariant ring of G_5_symp
S_quo, quo_map = quo(S, relations_ideal);

# Define the ideals corresponding to the polynomials that arise from the reflecting hyperplanes in the Bonnafé construction
delta_11_ideal = ideal(R, delta_11);
delta_12_ideal = ideal(R, delta_12);
delta_21_ideal = ideal(R, delta_21);
delta_22_ideal = ideal(R, delta_22);

# Investigate the first polynomial delta_11 and compute the variety along which the blow-up shall be performed
W_delta_11 = preimage(pi_1, delta_11_ideal);
basis_W_delta_11 = standard_basis(W_delta_11, ordering=negdegrevlex(S));
quo_W_delta_11 = quo_map(W_delta_11);
quo_basis_W_delta_11 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_11]);

nonzero_quo_W_delta_11 = gens(quo_basis_W_delta_11)[[1:5; 7:10; 12:14; 17]];
test_delta_11 = ideal(S_quo, nonzero_quo_W_delta_11[1:5]);
print(test_delta_11 == quo_W_delta_11);  # Returns true, therefore, the ideal image of W_delta_11 in the quotient ring is generated by the first 5 elements of the basis
rad_test_delta_11 = radical(test_delta_11);
print(rad_test_delta_11 == test_delta_11);  # Retuns true and thus shows that the ideal is radical, so blowing up along this ideal is equivalent to blowing up the corresponding subvariety

# List the generators of the ideal W_delta_11_prime which satisfies W_delta_11 = W_delta_11_prime + relations_ideal
w_11_1 = (2*zeta^4 - 1)*y[1]*y[2] + 3*y[5];
w_11_2 = 6*y[1]^4 + (8*zeta^4 - 4)*y[1]*y[3] + 3*y[6];
w_11_3 = (-2*zeta^4 + 1)*y[2]^2 + 3*y[8];
w_11_4 = (-4*zeta^4 + 2)*y[1]^2*y[5] + (-2*zeta^4 + 1)*y[2]*y[3] + 3*y[9];
w_11_5 = -64*y[1]^3*y[3] + (32//3*zeta^4 - 16//3)*y[2]*y[4] + (-128//3*zeta^4 + 64//3)*y[3]^2 + 48*y[10];

W_delta_11_prime = ideal(S, [w_11_1, w_11_2, w_11_3, w_11_4, w_11_5]);
print(W_delta_11 == W_delta_11_prime + relations_ideal);  # Returns true

# Investigate the second polynomial delta_12 and check that it yields the same radical ideal as delta_11, which is expected as delta_12 = delta_11^2
W_delta_12 = preimage(pi_1, delta_12_ideal);
basis_W_delta_12 = standard_basis(W_delta_12, ordering=negdegrevlex(S));
quo_W_delta_12 = quo_map(W_delta_12);
quo_basis_W_delta_12 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_12]);

nonzero_quo_W_delta_12 = gens(quo_basis_W_delta_12)[[1:2; 6; 9; 11:13; 15; 18:19; 22; 25; 29:34; 39:41]];
test_delta_12 = ideal(S_quo, nonzero_quo_W_delta_12);
print(test_delta_12 == quo_W_delta_12);  # Returns true, therefore, the ideal image of W_delta_12 in the quotient ring is generated by the given elements of the basis
rad_test_delta_12 = radical(test_delta_12);
print(rad_test_delta_12 == test_delta_11);  # Returns true

w_12_1 = (-2*zeta^4 + 1)*y[2]^2 + 3*y[8];
w_12_2 = y[1]^3*y[2] + (-6*zeta^4 + 3)*y[1]^2*y[5] + (-2*zeta^4 + 1)*y[2]*y[3] + 3*y[9];
w_12_3 = (2*zeta^4 - 1)*y[1]^2*y[8] + (4*zeta^4 - 2)*y[1]*y[2]*y[5] + 3*y[5]^2;
w_12_4 = -2700*y[1]^4*y[5] - 2250*y[1]^2*y[2]*y[3] + (2700*zeta^4 - 1350)*y[1]^2*y[9] + (2250*zeta^4 - 1125)*y[1]*y[2]*y[6] + 675*y[5]*y[6];
w_12_5 = 216*y[1]^8 - 396//7*y[1]^2*y[2]*y[4] + (-864*zeta^4 + 432)*y[1]^2*y[10] + (-288*zeta^4 + 144)*y[1]*y[2]*y[7] + (576*zeta^4 - 288)*y[1]*y[3]*y[6] - 540//7*y[5]*y[7] + 918//7*y[6]^2

W_delta_12_prime = ideal(S, [w_12_1, w_12_2, w_12_3, w_12_4, w_12_5]);
print(W_delta_12 == W_delta_12_prime + relations_ideal);  # Returns true

# Investigate the third polynomial delta_21 and compute the variety along which the blow-up shall be performed
W_delta_21 = preimage(pi_1, delta_21_ideal);
basis_W_delta_21 = standard_basis(W_delta_21, ordering=negdegrevlex(S));
quo_W_delta_21 = quo_map(W_delta_21);
quo_basis_W_delta_21 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_21]);

nonzero_quo_W_delta_21 = gens(quo_basis_W_delta_21)[[1:5; 7:10; 12:14; 17]];
test_delta_21 = ideal(S_quo, nonzero_quo_W_delta_21[1:5]);
print(test_delta_21 == quo_W_delta_21);  # Returns true, therefore, the ideal image of W_delta_21 in the quotient ring is generated by the first 5 elements of the basis
rad_test_delta_21 = radical(test_delta_21);
print(rad_test_delta_21 == test_delta_21);

# List the generators of the ideal W_delta_21_prime which satisfies W_delta_21 = W_delta_21_prime + relations_ideal
w_21_1 = (-2*zeta^4 + 1)*y[1]*y[2] + 3*y[5];
w_21_2 =  6*y[1]^4 + (-8*zeta^4 + 4)*y[1]*y[3] + 3*y[6];
w_21_3 = (2*zeta^4 - 1)*y[2]^2 + 3*y[8];
w_21_4 = (4*zeta^4 - 2)*y[1]^2*y[5] + (2*zeta^4 - 1)*y[2]*y[3] + 3*y[9];
w_21_5 = -64*y[1]^3*y[3] + (-32//3*zeta^4 + 16//3)*y[2]*y[4] + (128//3*zeta^4 - 64//3)*y[3]^2 + 48*y[10];

W_delta_21_prime = ideal(S, [w_21_1, w_21_2, w_21_3, w_21_4, w_21_5]);
print(W_delta_21 == W_delta_21_prime + relations_ideal);  # Returns true

# Investigate the second polynomial delta_22 and check that it yields the same radical ideal as delta_21, which is expected as delta_22 = delta_21^2
W_delta_22 = preimage(pi_1, delta_22_ideal);
basis_W_delta_22 = standard_basis(W_delta_22, ordering=negdegrevlex(S));
quo_W_delta_22 = quo_map(W_delta_22);
quo_basis_W_delta_22 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_22]);

# Investigate the polynomial delta_11*delta_21
W_delta_11_delta_21 = preimage(pi_1, ideal(R, delta_11*delta_21));
basis_W_delta_11_delta_21 = standard_basis(W_delta_11_delta_21, ordering=negdegrevlex(S));
quo_W_delta_11_delta_21 = quo_map(W_delta_11_delta_21);
quo_basis_W_delta_11_delta_21 = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_delta_11_delta_21]);

nonzero_quo_W_delta_11_delta_21 = gens(quo_basis_W_delta_11_delta_21)[[2; 5:6; 9; 11:17; 20:21; 23:24; 27; 30:32; 35; 38; 40; 42:43; 46:48; 53:55]];
test_delta_11_delta_21 = ideal(S_quo, nonzero_quo_W_delta_11_delta_21[[1:5; 8:11; 20]]);
print(test_delta_11_delta_21 == quo_W_delta_11_delta_21);
rad_test_delta_11_delta_21 = radical(test_delta_11_delta_21);
print(rad_test_delta_11_delta_21 == test_delta_11_delta_21);  
print("test done");

w_11_21_1 = y[1]*y[8] + y[2]*y[5];
w_11_21_2 = y[1]^2*y[2]^2 + 3*y[5]^2;
w_11_21_3 = -2*y[1]^4*y[2] + 12*y[1]*y[9] + 3*y[2]*y[6];
w_11_21_4 = 6*y[1]^4*y[5] + 4*y[1]^2*y[2]*y[3] + 3*y[5]*y[6];
w_11_21_5 = 12*y[1]^8 - 2*y[1]^2*y[2]*y[4] + 6*y[5]*y[7] - 3*y[6]^2;
w_11_21_6 = -y[1]*y[2]^3 + 3*y[5]*y[8];
w_11_21_7 = -y[1]*y[2]^2*y[3] - 9*y[5]*y[9] + 3*y[6]*y[8];
w_11_21_8 = -y[1]*y[2]^2*y[4] - 18*y[5]*y[10] + 3*y[7]*y[8];
w_11_21_9 = y[2]^4 + 3*y[8]^2;
w_11_21_10 = 162*y[1]^4*y[10] + 36*y[1]^3*y[2]*y[7] + 18*y[1]*y[2]*y[3]*y[4] + 54*y[5]*y[11] - 27*y[6]*y[10];

W_delta_11_delta_21_prime = ideal(S, [w_11_21_1, w_11_21_2, w_11_21_3, w_11_21_4, w_11_21_5, w_11_21_6, w_11_21_7, w_11_21_8, w_11_21_9, w_11_21_10]);
print(W_delta_11_delta_21 == W_delta_11_delta_21_prime + relations_ideal);  # Returns true

# Compare this to the ideal in which Berry computes the blowup. This ideal is corresponds to the reduced singular locus
R44 = G22*G22_conj;
R93 = F31*F40_conj*G22_conj;
R93prime = F40*F31_conj*G22;
R66 = F04*F40_conj*G22_conj;
R66prime = F40*F04_conj*G22;
R39 = F04*F31_conj*G22_conj;
R39prime = F13*F04_conj*G22;
R131 = F40*F40_conj*G51;
R113 = F04*F04_conj*G15;
R142 = F40*F31*F40_conj*F31_conj;
R214 = F13*F04*F13_conj*F04_conj;
R191 = F40*F31*F40_conj^3;
R119 = F13*F04*F04_conj^3;
R240 = F40^3*F40_conj^3;
R024 = F04^3*F04_conj^3;

Berry_ideal = ideal(R, [R44, R93, R93prime, R66, R66prime, R39, R39prime, R131, R113, R142, R214, R191, R119, R240, R024]);

# Compute the ideal in S which corresponds to the ideal given by Berry
W_Berry_ideal = preimage(pi_1, Berry_ideal);
basis_W_Berry_ideal = standard_basis(W_Berry_ideal, ordering=negdegrevlex(S));
quo_W_Berry_ideal = quo_map(W_Berry_ideal);
quo_basis_W_Berry_ideal = ideal(S_quo, [simplify(quo_map(b)) for b in basis_W_Berry_ideal]);
print(quo_basis_W_Berry_ideal);

# Check if any of the previously computed ideals is equal the ideal corresponding to Berry's ideal
print(quo_basis_W_Berry_ideal == quo_W_delta_11_prime + relations_ideal);  # Returns 
print(quo_basis_W_Berry_ideal == quo_W_delta_21_prime + relations_ideal);  # Returns
print(quo_basis_W_Berry_ideal == quo_W_delta_11_delta_21_prime + relations_ideal);  # Returns
