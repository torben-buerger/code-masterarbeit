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

# Define the genarators of the ideals W_delta_11, W_delta_21, and W_delta_11_21 in which we want to compute the blowup via the Rees algebra
w_11_1 = y[7];
w_11_2 = -y[2]^3 + y[10];
w_11_3 = -3*y[1]*y[5] + y[2]*y[4];
w_11_4 = 4*y[1]^3*y[4] - y[2]*y[6] - 3*y[3]*y[5] + 4*q3*y[4]^2;
# Define the vector used to construct the underlying map of the blowup along the divisor W_delta_11
blowup_vector_11 = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], t[1]*w_11_1, t[1]*w_11_2, t[1]*w_11_3, t[1]*w_11_4];

w_21_1 = y[1]^2*y[2] - 3*q3*y[5];
w_21_2 = y[1]*y[2]^2 - 6*q3*y[7];
w_21_3 = -y[2]^3 + 2*y[10];
w_21_4 = -4*q3*y[1]*y[4] + y[2]*y[3];
# Define the vector used to construct the underlying map of the blowup along the divisor W_delta_21
blowup_vector_21 = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], t[1]*w_21_1, t[1]*w_21_2, t[1]*w_21_3, t[1]*w_21_4];

w_11_21_1 = 6*q3*y[1]^3*y[7] - 3*q3*y[1]^2*y[2]*y[5] + q3*y[1]*y[2]^2*y[4] + 18*y[4]*y[7];
w_11_21_2 = q3*y[1]*y[2]^3 - q3*y[1]*y[10] + 9*y[2]*y[7];
w_11_21_3 = q3*y[1]^2*y[2]*y[7] + 9*y[5]*y[7];
w_11_21_4 = 36*y[1]^2*y[4]*y[5] + 3*q3*y[1]*y[2]*y[3]*y[5] - 12*y[1]*y[2]*y[4]^2 - q3*y[2]^2*y[3]*y[4];
w_11_21_5 = q3*y[1]*y[2]^2*y[7] + 18*y[7]^2;
w_11_21_6 = -y[2]^3*y[7] + 2*y[7]*y[10];
w_11_21_7 = y[2]^6 - 3*y[2]^3*y[10] + 2*y[10]^2;
w_11_21_8 = -2*y[1]^4*y[2]*y[6] + 8*q3*y[1]^2*y[5]*y[6] - q3*y[1]*y[2]^2*y[8] - (11//3)q3*y[1]*y[2]*y[4]*y[6] - y[1]*y[3]^2*y[7] + -9*q3*y[1]*y[3]*y[4]*y[5] - 8*y[1]*y[4]^3 + 3//4*y[2]^2*y[3]*y[6] + 9//4*y[2]*y[3]^2*y[5] - 26*y[7]*y[8];
# Define the vector used to construct the underlying map of the blowup along the divisor W_delta_11_21
blowup_vector_11_21 = [y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], t[1]*w_11_21_1, t[1]*w_11_21_2, t[1]*w_11_21_3, t[1]*w_11_21_4, t[1]*w_11_21_5, t[1]*w_11_21_6, t[1]*w_11_21_7, t[1]*w_11_21_8];

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

# Compute the blowup along W_delta_21
blowup_map_21 = hom(C, Bt, blowup_vector_21);
kernel_blowup_map_21 = preimage(blowup_map_21, relations_ideal);
basis_blowup_21 = standard_basis(kernel_blowup_map_21, ordering=negdegrevlex(C));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_21 = ideal(C, collect(basis_blowup_21));
# Elimination theory approach to get the ideal of the blowup variety in the relative projective space P^4_C
blowup_ideal_21 = ideal(Aux_large_C, union([embed_B_1(rel) for rel in gens(relations_ideal)], [b[i] - embed_B_1(blowup_vector_21[10+i])*t[1] for i in 1:4]));
Rees_rel_21 = eliminate(blowup_ideal_21, [t[1]]);  # Get a presentation of the Rees algebra which corresponds to the blowup in delta_21


# Define the ring D=B[b_1,...,b_8] used to compute the rees algebra of the blowup in W_delta_11_21, which can be embedded into this relative projective space as the ideal along which we blow up are generated by 8 elements
D, y, b = polynomial_ring(M, :y => (1:10), :b => (1:8));
Aux_large_D, y, t, b = polynomial_ring(M, :y => (1:10), :t => (1:1), :b => (1:8));
embed_B_2 = hom(Bt, Aux_large_D, union([y[i] for i in 1:10], [t[1]]));
embed_D = hom(D, Aux_large_D, union([y[i] for i in 1:10], [b[i] for i in 1:8]));

# Compute the blowup along W_delta_11_21
blowup_map_11_21 = hom(D, Bt, blowup_vector_11_21);
kernel_blowup_map_11_21 = preimage(blowup_map_11_21, relations_ideal);
basis_blowup_11_21 = standard_basis(kernel_blowup_map_11_21, ordering=negdegrevlex(D));  # Get a smaller generating set of the blowup ideal
basis_blowup_ideal_11_21 = ideal(D, collect(basis_blowup_11_21));
# Elimination theory approach to get the ideal of the blowup variety in the relative projective space P^8_C
blowup_ideal_11_21 = ideal(Aux_large_D, union([embed_B_2(rel) for rel in gens(relations_ideal)], [b[i] - embed_B_2(blowup_vector_11_21[10+i])*t[1] for i in 1:8]));
Rees_rel_11_21 = eliminate(blowup_ideal_11_21, [t[1]]);  # Get a presentation of the Rees algebra which corresponds to the blowup in delta_11_21
