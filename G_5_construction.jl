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

conj = conjugacy_classes(G_5);

# Compute the corresponding symplectic reflection group
a_symp = block_diagonal_matrix([a, transpose(inv(a))]);
b_symp = block_diagonal_matrix([b, transpose(inv(b))]);
c_symp = block_diagonal_matrix([c, transpose(inv(c))]);
d_symp = block_diagonal_matrix([d, transpose(inv(d))]);
c_conj_symp = block_diagonal_matrix([c_conj, transpose(inv(c_conj))]);
G_5_symp = matrix_group([a_symp, b_symp, c_symp, c_conj_symp]);

# Compute the invariants of the symplectic reflection group
IR = invariant_ring(G_5_symp);
poly_invar = polynomial_ring(IR);
u = gens(poly_invar);
invars = fundamental_invariants(IR);
invars_ideal_1 = ideal(poly_invar, invars);

# In order to compare with the invariants given by Berry, we consider the corresponding polynomials in the polynomial ring R in 4 variables
R, x = polynomial_ring(K, :x => (1:4));

# Invariants of G_5_symp by the above computation
invar_1 = x[1]*x[3] + x[2]*x[4];
invar_2 = x[1]^5*x[2] - x[1]*x[2]^5;
invar_3 = x[1]*x[2]^2*x[3]^3 - x[2]^3*x[3]^2*x[4] - x[1]^3*x[3]*x[4]^2 + x[1]^2*x[2]*x[4]^3;
invar_4 = x[3]^5*x[4] - x[3]*x[4]^5;
invar_5 = x[1]^4*x[2]^3*x[3] + 1//7*x[2]^7*x[3] - 1//7*x[1]^7*x[4] - x[1]^3*x[2]^4*x[4];
invar_6 = x[2]^4*x[3]^4 - 4*x[1]^3*x[2]*x[3]^3*x[4] + 6*x[1]^2*x[2]^2*x[3]^2*x[4]^2 - 4*x[1]*x[2]^3*x[3]*x[4]^3 + x[1]^4*x[4]^4;
invar_7 = x[2]*x[3]^7 - 7*x[1]*x[3]^4*x[4]^3 + 7*x[2]*x[3]^3*x[4]^4 - x[1]*x[4]^7;
invar_8 = x[3]^12 - 33*x[3]^8*x[4]^4 - 33*x[3]^4*x[4]^8 + x[4]^12;
invar_9 = x[2]^3*x[3]^9 + 3//5*x[1]^3*x[3]^8*x[4] - 36//5*x[1]^2*x[2]*x[3]^7*x[4]^2 + 84//5*x[1]*x[2]^2*x[3]^6*x[4]^3 - 42//5*x[2]^3*x[3]^5*x[4]^4 + 42//5*x[1]^3*x[3]^4*x[4]^5 - 84//5*x[1]^2*x[2]*x[3]^3*x[4]^6 + 36//5*x[1]*x[2]^2*x[3]^2*x[4]^7 -3//5*x[2]^3*x[3]*x[4]^8 - x[1]^3*x[4]^9;
invar_10 = x[1]^4*x[2]^2*x[3]^6 - 1//17*x[2]^6*x[3]^6 + 24//17*x[1]^3*x[2]^3*x[3]^5*x[4] - 15//17*x[1]^6*x[3]^4*x[4]^2 + 15//17*x[1]^2*x[2]^4*x[3]^4*x[4]^2 - 40//17*x[1]^5*x[2]*x[3]^3*x[4]^3 - 40//17*x[1]*x[2]^5*x[3]^3*x[4]^3 + 15//17*x[1]^4*x[2]^2*x[3]^2*x[4]^4 - 15//17*x[2]^6*x[3]^2*x[4]^4 + 24//17*x[1]^3*x[2]^3*x[3]*x[4]^5 - 1//17*x[1]^6*x[4]^6 + x[1]^2*x[2]^4*x[4]^6;
invar_11 = x[1]^8*x[2]*x[3]^3 - 14//23*x[1]^4*x[2]^5*x[3]^3 - 1//23*x[2]^9*x[3]^3 + 60//23*x[1]^7*x[2]^2*x[3]^2*x[4] - 84//23*x[1]^3*x[2]^6*x[3]^2*x[4] + 84//23*x[1]^6*x[2]^3*x[3]*x[4]^2 - 60//23*x[1]^2*x[2]^7*x[3]*x[4]^2 + 1//23*x[1]^9*x[4]^3 + 14//23*x[1]^5*x[2]^4*x[4]^3 - x[1]*x[2]^8*x[4]^3;
invar_12 = x[1]^12 - 33*x[1]^8*x[2]^4 - 33*x[1]^4*x[2]^8 + x[2]^12;
invars_ideal_1 = ideal(R, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, invar_9, invar_10, invar_11, invar_12]);

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

gen_ideal = ideal(R, [F11, F40, F31, G22, F13, F04, F60, G51, G42, F33, G24, G15, F06, F40_conj, F31_conj, G22_conj, F13_conj, F04_conj]);

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

invars_ideal_2 = ideal(R, [h, A1, A0, Aneg1, B1, B0, Bneg1, C2, C1, C0, Cneg1, Cneg2]);
print(invars_ideal_1 == invars_ideal_2)  # confirm that both ideals are equal

# Now express the invariants invars_ideal_1 in terms of the invariants and eigenvectors given by Berry
I_9 = ideal(R, [Aneg1, Cneg1]);
print(ideal_membership(invar_9, I_9));
coeffs_9 = coordinates(invar_9, I_9);
print(1//10*F11^3*Aneg1+1//2*Cneg1 == invar_9);

I_10 = ideal(R, [A0, C0]);
print(ideal_membership(invar_10, I_10));
coeffs_10 = coordinates(invar_10, I_10);
print(-16//51*F11^3*A0-1//34*C0 == invar_10);

I_11 = ideal(R, [A1, C1]);
print(ideal_membership(invar_11, I_11));
coeffs_11 = coordinates(invar_11, I_11);
print(1//6*F11^3*A1-1//46*C1 == invar_11);

# The rest of the invariants can be directly expressed as scalar multiples, thus we can work with the invariants given by Berry

# Define polynomial ring in 8 variables to express the relations between the invariants
S, y = polynomial_ring(K, :y => (1:12));
# Define homomorphism from S to R mapping y_i to the invariants given by Berry
f = hom(S, R, [h, A1, A0, Aneg1, B1, B0, Bneg1, C2, C1, C0, Cneg1, Cneg2]);
# Compute kernel of f to obtain the relations
relations_ideal = kernel(f);
basis_relations = standard_basis(relations_ideal, ordering=negdegrevlex(S));

g = hom(S, R, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, invar_9, invar_10, invar_11, invar_12]);
relations_ideal_1 = kernel(g);
print(relations_ideal == relations_ideal_1)  # confirm that both ideals are equal