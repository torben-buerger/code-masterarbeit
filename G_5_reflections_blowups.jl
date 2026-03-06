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

# Find the reflecting hyperplanes (= lines as we are in dimension 2) of the reflections
ref_1_eigvecs = eigenspaces(K, matrix(ref_1));
ref_2_eigvecs = eigenspaces(K, matrix(ref_2));
ref_3_eigvecs = eigenspaces(K, matrix(ref_3));
ref_4_eigvecs = eigenspaces(K, matrix(ref_4));
