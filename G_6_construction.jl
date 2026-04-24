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

S, y = polynomial_ring(K, :y => (1:10));
pi = hom(S, R, [gens(invars_ideal_2)[i] for i in 1:10]);
relations_ideal = kernel(pi);
relations_basis = standard_basis(relations_ideal, ordering=negdegrevlex(S));
print(length(relations_basis));
