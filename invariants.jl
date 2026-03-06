using Oscar
Qx, s = QQ["s"];

# For the rest of the comptations, we only need the algebraic number sqrt(-3), hence we can pass to a smaller number field
M, q3 = number_field(s^2 + 3, "q3");
A, x = polynomial_ring(M, :x => (1:4));

# Define the invariant polynomials in this smaller number field, see matrixgroups.jl for the derivation, we have to be careful with the order to be able to usee the same simplified relations as in the paper of Lehn and Sorger
invar_1 = x[1]*x[3] + x[2]*x[4];
invar_2 = x[3]^4 - 2*q3*x[3]^2*x[4]^2 + x[4]^4;
invar_3 = x[1]^4 + 2*q3*x[1]^2*x[2]^2 + x[2]^4;
invar_4 = x[2]*x[3]^3 - q3*x[1]*x[3]^2*x[4] + q3*x[2]*x[3]*x[4]^2 - x[1]*x[4]^3;
invar_5 = -q3*x[1]^2*x[2]*x[3] + x[2]^3*x[3] - x[1]^3*x[4] + q3*x[1]*x[2]^2*x[4];
invar_6 = x[1]^5*x[2] - x[1]*x[2]^5;
invar_7 = x[3]^5*x[4] - x[3]*x[4]^5;
invar_8 = x[1]*x[2]^2*x[3]^3 - x[2]^3*x[3]^2*x[4] - x[1]^3*x[3]*x[4]^2 + x[1]^2*x[2]*x[4]^3;
invars_ideal = ideal(A, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8]);

# Define polynomial ring in 8 variables to express the relations between the invariants
B, z = polynomial_ring(M, :z => (1:8));
# Define homomorphism from B to A mapping z_i to invar_i
f = hom(B, A, [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8]);
# Compute kernel of f to obtain the relations
relations_ideal = kernel(f);
basis_relations = standard_basis(relations_ideal, ordering=negdegrevlex(B));
#print(basis_relations)

# Define relations and reorder/normalize them if possible to match the relations given in the paper of Lehn and Sorger
rel_1 = q3*z[1]^3*z[5] - z[1]*z[3]*z[4] - 2*z[2]*z[6] - z[5]*z[8];  # == -basis_relations[2]
rel_2 = q3*z[1]^3*z[4] + z[1]*z[2]*z[5] - 2*z[3]*z[7] - z[4]*z[8];  # == -basis_relations[4]
rel_3 = z[1]*z[5]^2 + z[3]*z[8] + 2*z[4]*z[6];  # == basis_relations[3]
rel_4 = z[1]*z[4]^2 - z[2]*z[8] - 2*z[5]*z[7];  # == -basis_relations[5]
rel_5 = -z[1]^4 - 3*q3*z[1]*z[8] + z[2]*z[3] - z[4]*z[5]; # == basis_relations[1]
rel_6 = z[1]^3*z[8] + 1//4*z[2]*z[5]^2 - 1//4*z[3]*z[4]^2 - q3*z[6]*z[7] + q3*z[8]^2; # == -1//12*q3*basis_relations[6]
rel_7 = -2*z[1]^3*z[6] + q3*z[1]^2*z[3]*z[5] - z[3]^2*z[4] + z[5]^3 - 6*q3*z[6]*z[8];  # == -1//3*q3*basis_relations[7]
rel_8 = -2*z[1]^3*z[7] + q3*z[1]^2*z[2]*z[4] + z[2]^2*z[5] - z[4]^3 - 6*q3*z[7]*z[8];  # == -1//3*q3*basis_relations[8]
rel_9 = z[1]^3*z[8] + q3*z[1]^2*z[4]*z[5] + z[2]*z[5]^2 - z[3]*z[4]^2 + 3*q3*z[8]^2;  # == 1//3*q3*basis_relations[9]
relations_ideal_1 = ideal(B, [rel_1, rel_2, rel_3, rel_4, rel_5, rel_6, rel_7, rel_8, rel_9]);
print(relations_ideal == relations_ideal_1)  # confirm that both ideals are equal

# Define the relations as given in the paper of Lehn and Sorger, which differ slightly from the computed ones
rel_paper_1 = q3*z[1]^3*z[5] - z[1]*z[3]*z[4] - 2*z[2]*z[6] - z[5]*z[8];
rel_paper_2 = q3*z[1]^3*z[4] + z[1]*z[2]*z[5] - 2*z[3]*z[7] - z[4]*z[8];
rel_paper_3 = z[1]*z[5]^2 + 2*z[4]*z[6] + z[3]*z[8];
rel_paper_4 = z[1]*z[4]^2 - 2*z[5]*z[7] - z[2]*z[8];
rel_paper_5 = -z[1]^4 + z[2]*z[3] - z[4]*z[5] - 3*q3*z[1]*z[8];
rel_paper_6 = z[1]^2*z[4]*z[5] + q3*z[1]^3*z[8] + 4*z[6]*z[7] - z[8]^2;
rel_paper_7 = q3*z[1]^2*z[3]*z[5] - 2*z[1]^3*z[6] - z[3]^2*z[4] + z[5]^3 - 6*q3*z[6]*z[8];
rel_paper_8 = q3*z[1]^2*z[2]*z[4] - 2*z[1]^3*z[7] - z[4]^3 + z[2]^2*z[5] - 6*q3*z[7]*z[8];
rel_paper_9 = 4*z[1]^2*z[4]*z[5] + q3*z[3]*z[4]^2 - q3*z[2]*z[5]^2 + 4*z[6]*z[7] + 8*z[8]^2;
relations_ideal_2 = ideal(B, [rel_paper_1, rel_paper_2, rel_paper_3, rel_paper_4, rel_paper_5, rel_paper_6, rel_paper_7, rel_paper_8, rel_paper_9]);
print(relations_ideal == relations_ideal_2)  # confirm that both ideals are equal

#=
# Define the invariants of the dual representation in the smaller number field
invar_1_dual = x[1]*x[3] + x[2]*x[4];
invar_2_dual = x[1]^4 + (-2*q3)*x[1]^2*x[2]^2 + x[2]^4;
invar_3_dual = q3*x[1]^2*x[2]*x[3] + x[2]^3*x[3] - x[1]^3*x[4] + (-q3)*x[1]*x[2]^2*x[4];
invar_4_dual = x[2]*x[3]^3 + q3*x[1]*x[3]^2*x[4] + (-q3)*x[2]*x[3]*x[4]^2 - x[1]*x[4]^3;
invar_5_dual = x[3]^4 + (2*q3)*x[3]^2*x[4]^2 + x[4]^4;
invar_6_dual = x[1]^5*x[2] - x[1]*x[2]^5;
invar_7_dual = x[1]*x[2]^2*x[3]^3 - x[2]^3*x[3]^2*x[4] - x[1]^3*x[3]*x[4]^2 + x[1]^2*x[2]*x[4]^3;
invar_8_dual = x[3]^5*x[4] - x[3]*x[4]^5;

invars_dual_ideal = ideal(A, [invar_1_dual, invar_2_dual, invar_3_dual, invar_4_dual, invar_5_dual, invar_6_dual, invar_7_dual, invar_8_dual]);
# Define homomorphism from B to A mapping z_i to invar_i_dual
f_dual = hom(B, A, [invar_1_dual, invar_2_dual, invar_3_dual, invar_4_dual, invar_5_dual, invar_6_dual, invar_7_dual, invar_8_dual]);
# Compute kernel of f_dual to obtain the relations
relations_dual_ideal = kernel(f_dual);
=#

# Define the polynomials definining the vanishing locus of the union of the 4 lines given by the eigenvectors in the smaller number field (see matrixgroups.jl for the derivation)
C_1 = x[1]^4 - 2*q3*x[1]^2*x[2]^2 + x[2]^4;
C_2 = x[3]^4 + 2*q3*x[3]^2*x[4]^2 + x[4]^4;

C_1_ideal = ideal(A, C_1);
C_2_ideal = ideal(A, C_2);

# Determine the generators of the Weil divisors in the quotient ring B / relations_ideal_1 which is isomorphic to the invariant ring
B_quo, pi = quo(B, relations_ideal_1);

# Compute Weil divisor W_1 by taking preimages of the ideals under the homomorphism f
W_1 = preimage(f, C_1_ideal);
basis_W_1 = standard_basis(W_1, ordering=negdegrevlex(B));

#print([simplify(pi(f)) for f in basis_W_1]) # Per inspection of the basis elements, we can drop the elements 6, 9, 10, 14 from the basis since they are 0 in the quotient ring
# We also see that the basis elements 3 and 12 are redundant in the quotient ring, because they are scalar multiples of the 1. resp. the 5. basis element
#test_quo_ideal_1 = ideal(B_quo, [pi(basis_W_1[1]), pi(basis_W_1[2]), pi(basis_W_1[4])]);
#ideal_membership(pi(basis_W_1[11]), test_quo_ideal_1)  # Returns true for 11, 13, 15, 16, 17, hence we can drop these elements as well

# Define the generators of W_1 in the quotient ring B / relations_ideal_1, normalized and reordered to match the ones for W_2
W_1_gen_1 = z[2]*z[6] + 2*z[5]*z[8];  # == basis_W_1[5]
W_1_gen_2 = z[3]*z[5] + 2*q3*z[1]*z[6];  # == basis_W_1[2]
W_1_gen_3 = z[2]*z[3] - 4*q3*z[1]*z[8];  # == basis_W_1[1]
W_1_gen_4 = z[3]^3 - 12*q3*z[6]^2;  # == 1//3*q3*basis_W_1[8]
W_1_gen_5 = z[1]*z[3]^2 + 6*z[5]*z[6];  # == -basis_W_1[8]
W_1_gen_6 = z[1]^2*z[3] + q3*z[5]^2;  # == -1//q3*basis_W_1[3]

W_1_ideal = ideal(B, [W_1_gen_1, W_1_gen_2, W_1_gen_3, W_1_gen_4, W_1_gen_5, W_1_gen_6]) + relations_ideal_1;
print(W_1 == W_1_ideal)  # confirm that both ideals are equal

# As we want to blow up the quotient variety in the variety defined by W_2, we need to confirm that the ideal of W_2 in the quotient ring is radical
W_1_quo_ideal = ideal(B_quo, [pi(W_1_gen_1), pi(W_1_gen_2), pi(W_1_gen_3), pi(W_1_gen_4), pi(W_1_gen_5), pi(W_1_gen_6)]);
print(radical(W_1_quo_ideal) == W_1_quo_ideal)  # confirm that the ideal is radical, hence it is the vanishing ideal of the Weil divisor W_2


# Compute Weil divisor W_2 by taking preimages of the ideals under the homomorphism f
W_2 = preimage(f, C_2_ideal);
basis_W_2 = standard_basis(W_2, ordering=negdegrevlex(B));

#print([simplify(pi(f)) for f in basis_W_2]) # Per inspection of the basis elements, we can drop the elements 5, 6, 9, 13 from the basis since they are 0 in the quotient ring
# We also see that the basis elements 4 and 12 are redundant in the quotient ring, because they are scalar multiples of the 1. resp. the 7. basis element
#test_quo_ideal_2 = ideal(B_quo, [pi(basis_W_2[1]), pi(basis_W_2[2]), pi(basis_W_2[3])]);
#ideal_membership(pi(basis_W_2[10]), test_quo_ideal_2)  # Returns true for 10, 14, 15, 16, 17, hence we can drop these elements as well

# Define the generators of W_2 in the quotient ring B / relations_ideal_1, normalized and reordered to match the ones given in the paper of Lehn and Sorger
W_2_gen_1 = z[3]*z[7] + 2*z[4]*z[8];  # == basis_W_2[7]
W_2_gen_2 = z[2]*z[4] + 2*q3*z[1]*z[7];  # == basis_W_2[2]
W_2_gen_3 = z[2]*z[3] - 4*q3*z[1]*z[8];  # == basis_W_2[1]
W_2_gen_4 = z[2]^3 + 12*q3*z[7]^2;  # == 1//3*q3*basis_W_2[11]
W_2_gen_5 = z[1]*z[2]^2 - 6*z[4]*z[7];  # == -basis_W_2[8]
W_2_gen_6 = z[1]^2*z[2] - q3*z[4]^2;  # == 1//q3*basis_W_2[3]

W_2_ideal = ideal(B, [W_2_gen_1, W_2_gen_2, W_2_gen_3, W_2_gen_4, W_2_gen_5, W_2_gen_6]) + relations_ideal_1;
print(W_2 == W_2_ideal)  # confirm that both ideals are equal

# As we want to blow up the quotient variety in the variety defined by W_2, we need to confirm that the ideal of W_2 in the quotient ring is radical
W_2_quo_ideal = ideal(B_quo, [pi(W_2_gen_1), pi(W_2_gen_2), pi(W_2_gen_3), pi(W_2_gen_4), pi(W_2_gen_5), pi(W_2_gen_6)]);
print(radical(W_2_quo_ideal) == W_2_quo_ideal)  # confirm that the ideal is radical, hence it is the vanishing ideal of the Weil divisor W_2