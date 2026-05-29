// Define the base field and number field
Q := RationalField();
R_t<t> := PolynomialRing(Q);
K<i> := NumberField(t^2 + 1);      // Adjoin i = sqrt(-1)
K_z<z> := PolynomialRing(K);
L<q2> := NumberField(z^2 - 2);     // Adjoin q2 = sqrt(2)
i := L ! i;                        // Coerce i into the top field L

// Define the matrix group G_12
r3 := (1/q2) * Matrix(L, 2, 2, [1, -1, -1, -1]);
r3prime := (1/q2) * Matrix(L, 2, 2, [1, 1, 1, -1]);
r3primeprime := (1/q2) * Matrix(L, 2, 2, [0, 1+i, 1-i, 0]);
G_12 := MatrixGroup<2, L | r3, r3prime, r3primeprime>;

// Define the symplectic reflection group G_12_symp
r3_symp := DiagonalJoin(r3, Transpose(r3^-1));
r3prime_symp := DiagonalJoin(r3prime, Transpose(r3prime^-1));
r3primeprime_symp := DiagonalJoin(r3primeprime, Transpose(r3primeprime^-1));
G_12_symp := MatrixGroup<4, L | r3_symp, r3prime_symp, r3primeprime_symp>;

// Compute the invariants of the symplectic reflection group automatically
IR := InvariantRing(G_12_symp);
invars_auto := FundamentalInvariants(IR);
invars_ideal_1 := ideal<PolynomialRing(IR) | invars_auto>;

// Define the manual invariants as polynomials in the original variables
R<x1, x2, x3, x4> := PolynomialRing(L, 4);
x := [x1, x2, x3, x4];

invar_1 := x1*x3 + x2*x4;
invar_2 := x1^5*x2 - x1*x2^5;
invar_3 := x1^4*x3*x4 - 2*x1^3*x2*x4^2 + 2*x1*x2^3*x3^2 - x2^4*x3*x4;
invar_4 := x1^2*x3*x4^3 + 1/2*x1*x2*x3^4 - 1/2*x1*x2*x4^4 - x2^2*x3^3*x4;
invar_5 := x3^5*x4 - x3*x4^5;
invar_6 := x1^8 + 14*x1^4*x2^4 + x2^8;
invar_7 := x1^6*x4^2 + 3*x1^4*x2^2*x3^2 - 8*x1^3*x2^3*x3*x4 + 3*x1^2*x2^4*x4^2 + x2^6*x3^2;
invar_8 := x1^4*x4^4 - 4*x1^3*x2*x3^3*x4 + 6*x1^2*x2^2*x3^2*x4^2 - 4*x1*x2^3*x3*x4^3 + x2^4*x3^4;
invar_9 := x1^2*x3^4*x4^2 + 1/3*x1^2*x4^6 - 8/3*x1*x2*x3^3*x4^3 + 1/3*x2^2*x3^6 + x2^2*x3^2*x4^4;
invar_10 := x3^8 + 14*x3^4*x4^4 + x4^8;
invar_11 := x1^11*x4 + 11*x1^8*x2^3*x3 - 22*x1^7*x2^4*x4 + 22*x1^4*x2^7*x3 - 11*x1^3*x2^8*x4 - x2^11*x3;
invar_12 := x1^9*x4^3 - 9*x1^7*x2^2*x3^2*x4 + 15*x1^6*x2^3*x3*x4^2 - 9*x1^5*x2^4*x4^3 + 9*x1^4*x2^5*x3^3 - 15*x1^3*x2^6*x3^2*x4 +
            9*x1^2*x2^7*x3*x4^2 - x2^9*x3^3;
invar_13 := x1^7*x4^5 + 5*x1^6*x2*x3^3*x4^2 - 15*x1^5*x2^2*x3^2*x4^3 + 3*x1^4*x2^3*x3^5 + 10*x1^4*x2^3*x3*x4^4 - 10*x1^3*x2^4*x3^4*x4
        - 3*x1^3*x2^4*x4^5 + 15*x1^2*x2^5*x3^3*x4^2 - 5*x1*x2^6*x3^2*x4^3 - x2^7*x3^5;
invar_14 := x1^5*x4^7 + 3/2*x1^4*x2*x3^7 + 35/2*x1^4*x2*x3^3*x4^4 - 7/2*x1^3*x2^2*x3^6*x4 - 21/2*x1^3*x2^2*x3^2*x4^5 +
        21/2*x1^2*x2^3*x3^5*x4^2 + 7/2*x1^2*x2^3*x3*x4^6 - 35/2*x1*x2^4*x3^4*x4^3 - 3/2*x1*x2^4*x4^7 - x2^5*x3^7;
invar_15 := x1^3*x3^8*x4 - 1/9*x1^3*x4^9 + 2*x1^2*x2*x3^7*x4^2 - 14/3*x1^2*x2*x3^3*x4^6 + 14/3*x1*x2^2*x3^6*x4^3 - 2*x1*x2^2*x3^2*x4^7
        + 1/9*x2^3*x3^9 - x2^3*x3*x4^8;
invar_16 := x1*x3^8*x4^3 + 2*x1*x3^4*x4^7 - 1/11*x1*x4^11 + 1/11*x2*x3^11 - 2*x2*x3^7*x4^4 - x2*x3^3*x4^8;

invars_list := [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, 
               invar_9, invar_10, invar_11, invar_12, invar_13, invar_14, invar_15, invar_16];
invars_ideal_2 := ideal<R | invars_list>;

// Compute the relations between the invariants to obtain the presentation using elimination theory.
// Define a big polynomial ring where the first 4 variables will be eliminated
T<tx1,tx2,tx3,tx4, ty1,ty2,ty3,ty4,ty5,ty6,ty7,ty8,ty9,ty10,ty11,ty12,ty13,ty14,ty15,ty16> := PolynomialRing(L, 20, "elim", 4);

// Map the invariants from R into the x-variables of T
map_R_to_T := hom< R -> T | [tx1, tx2, tx3, tx4] >;
invars_in_T := [map_R_to_T(f) : f in invars_list];
ty_vars := [ty1, ty2, ty3, ty4, ty5, ty6, ty7, ty8, ty9, ty10, ty11, ty12, ty13, ty14, ty15, ty16];

// Build the graph ideal: (ty_i - invar_i) using correct bracketed range [1..16]
graph_ideal := ideal< T | [ty_vars[i] - invars_in_T[i] : i in [1..16]] >;

// Eliminate tx1, tx2, tx3, tx4
elim_ideal := EliminationIdeal(graph_ideal, 4);

// Define the target ring S with the desired "grevlexw" ordering working like the ordering available in Oscar
W := [Degree(f) : f in invars_list];
S<y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16> := PolynomialRing(L, 16, "grevlexw", W);

// Map the eliminated relations from T into S (sending x-variables to 0)
map_T_to_S := hom< T -> S | [0,0,0,0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16] >;
relations_ideal := ideal< S | [map_T_to_S(f) : f in Basis(elim_ideal)] >;

// Extract the relations basis under the specific ordering
relations_basis := Basis(relations_ideal);
print "Number of relations found:", #relations_basis;

relations_basis_ideal := ideal<S | relations_basis>;