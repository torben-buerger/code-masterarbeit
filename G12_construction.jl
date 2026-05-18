using Oscar
# In order to define the generators of G12 provided by Lehrer--Taylor, we need the algebraic numbers i = sqrt(-1) and q2 = sqrt(2)
K, i = cyclotomic_field(4, "i");
Kt, t = polynomial_ring(K, :t);
L, q2 = number_field(t^2-2, "q2");
i = L(i);
    
r3 = 1//q2 * matrix(L, 2, 2, [1 -1; -1 -1]);
r3prime = 1//q2 * matrix(L, 2, 2, [1 1; 1 -1]);
r3primeprime = 1//q2 * matrix(L, 2, 2, [0 1+i; 1-i 0]);

G_12 = matrix_group([r3, r3prime, r3primeprime]);
describe(G_12);

r3_symp = block_diagonal_matrix([r3, transpose(inv(r3))]);
r3prime_symp = block_diagonal_matrix([r3prime, transpose(inv(r3prime))]);
r3primeprime_symp = block_diagonal_matrix([r3primeprime, transpose(inv(r3primeprime))]);
G_12_symp = matrix_group([r3_symp, r3prime_symp, r3primeprime_symp]);

# Compute the invariants of the symplectic reflection group
IR = invariant_ring(G_12_symp);
poly_invar = polynomial_ring(IR);
u = gens(poly_invar);
invars = fundamental_invariants(IR);
invars_ideal_1 = ideal(poly_invar, invars);

# Define the invariants as polynomials in the original variables
R, x = polynomial_ring(K, :x => (1:4));

invar_1 = x[1]*x[3] + x[2]*x[4];
invar_2 = x[1]^5*x[2] - x[1]*x[2]^5;
invar_3 = x[1]*x[2]^3*x[3]^2 + 1//2*x[1]^4*x[3]*x[4] - 1//2*x[2]^4*x[3]*x[4] - x[1]^3*x[2]*x[4]^2;
invar_4 = x[1]*x[2]*x[3]^4 - 2*x[2]^2*x[3]^3*x[4] + 2*x[1]^2*x[3]*x[4]^3 - x[1]*x[2]*x[4]^4;
invar_5 = x[3]^5*x[4] - x[3]*x[4]^5;
invar_6 = x[1]^8 + 14*x[1]^4*x[2]^4 + x[2]^8;
invar_7 = x[1]^4*x[2]^2*x[3]^2 + 1//3*x[2]^6*x[3]^2 - 8//3*x[1]^3*x[2]^3*x[3]*x[4] + 1//3*x[1]^6*x[4]^2 + x[1]^2*x[2]^4*x[4]^2;
invar_8 = x[2]^4*x[3]^4 - 4*x[1]^3*x[2]*x[3]^3*x[4] + 6*x[1]^2*x[2]^2*x[3]^2*x[4]^2 - 4*x[1]*x[2]^3*x[3]*x[4]^3 + x[1]^4*x[4]^4;
invar_9 = x[2]^2*x[3]^6 + 3*x[1]^2*x[3]^4*x[4]^2 - 8*x[1]*x[2]*x[3]^3*x[4]^3 + 3*x[2]^2*x[3]^2*x[4]^4 + x[1]^2*x[4]^6;
invar_10 = x[3]^8 + 14*x[3]^4*x[4]^4 + x[4]^8;
invar_11 = x[2]*x[3]^11 + 11*x[1]*x[3]^8*x[4]^3 - 22*x[2]*x[3]^7*x[4]^4 + 22*x[1]*x[3]^4*x[4]^7 - 11*x[2]*x[3]^3*x[4]^8 - x[1]*x[4]^11;
invar_12 = x[2]^3*x[3]^9 + 3//5*x[1]^3*x[3]^8*x[4] - 36//5*x[1]^2*x[2]*x[3]^7*x[4]^2 + 84//5*x[1]*x[2]^2*x[3]^6*x[4]^3 - 42//5*x[2]^3*x[3]^5*x[4]^4 + 42//5*x[1]^3*x[3]^4*x[4]^5 - 84//5*x[1]^2*x[2]*x[3]^3*x[4]^6 + 36//5*x[1]*x[2]^2*x[3]^2*x[4]^7 - 3//5*x[2]^3*x[3]*x[4]^8 - x[1]^3*x[4]^9;
invar_13 = x[1]^4*x[2]*x[3]^7 - 1//19*x[2]^5*x[3]^7 + 14//19*x[1]^3*x[2]^2*x[3]^6*x[4] - 42//19*x[1]^2*x[2]^3*x[3]^5*x[4]^2 + 35//19*x[1]^5*x[3]^4*x[4]^3 - 105//19*x[1]*x[2]^4*x[3]^4*x[4]^3 + 105//19*x[1]^4*x[2]*x[3]^3*x[4]^4 - 35//19*x[2]^5*x[3]^3*x[4]^4 + 42//19*x[1]^3*x[2]^2*x[3]^2*x[4]^5 - 14//19*x[1]^2*x[2]^3*x[3]*x[4]^6 + 1//19*x[1]^5*x[4]^7 - x[1]*x[2]^4*x[4]^7;
invar_14 = x[1]^4*x[2]^3*x[3]^5 - 1//13*x[2]^7*x[3]^5 + 5//13*x[1]^7*x[3]^4*x[4] + 15//13*x[1]^3*x[2]^4*x[3]^4*x[4] + 10//13*x[1]^6*x[2]*x[3]^3*x[4]^2 + 30//13*x[1]^2*x[2]^5*x[3]^3*x[4]^2 - 30//13*x[1]^5*x[2]^2*x[3]^2*x[4]^3 - 10//13*x[1]*x[2]^6*x[3]^2*x[4]^3 - 15//13*x[1]^4*x[2]^3*x[3]*x[4]^4 - 5//13*x[2]^7*x[3]*x[4]^4 + 1//13*x[1]^7*x[4]^5 - x[1]^3*x[2]^4*x[4]^5;
invar_15 = x[1]^8*x[2]*x[3]^3 - 14//23*x[1]^4*x[2]^5*x[3]^3 - 1//23*x[2]^9*x[3]^3 + 60//23*x[1]^7*x[2]^2*x[3]^2*x[4] - 84//23*x[1]^3*x[2]^6*x[3]^2*x[4] + 84//23*x[1]^6*x[2]^3*x[3]*x[4]^2 - 60//23*x[1]^2*x[2]^7*x[3]*x[4]^2 + 1//23*x[1]^9*x[4]^3 + 14//23*x[1]^5*x[2]^4*x[4]^3 - x[1]*x[2]^8*x[4]^3;
invar_16 = x[1]^8*x[2]^3*x[3] + 2*x[1]^4*x[2]^7*x[3] - 1//11*x[2]^11*x[3] + 1//11*x[1]^11*x[4] - 2*x[1]^7*x[2]^4*x[4] - x[1]^3*x[2]^8*x[4];
invars_ideal = ideal(R, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, invar_9, invar_10, invar_11, invar_12, invar_13, invar_14, invar_15, invar_16]);

# Compute the relations between the invariants to obtain the presentation of the invariant ring
S, y = polynomial_ring(K, :y => (1:16));
pi = hom(S, R, [gens(invars_ideal)[i] for i in 1:16]);
relations_ideal = kernel(pi);
relations_basis = standard_basis(relations_ideal, ordering=negdegrevlex(S));
print(length(relations_basis));
relations_basis_ideal = ideal(S, relations_basis);
