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

classes = conjugacy_classes(G_5);

# Compute the corresponding symplectic reflection group
a_symp = block_diagonal_matrix([a, transpose(inv(a))]);
b_symp = block_diagonal_matrix([b, transpose(inv(b))]);
c_symp = block_diagonal_matrix([c, transpose(inv(c))]);
d_symp = block_diagonal_matrix([d, transpose(inv(d))]);
c_conj_symp = block_diagonal_matrix([c_conj, transpose(inv(c_conj))]);
G_5_symp = matrix_group([a_symp, b_symp, c_symp, c_conj_symp]);
