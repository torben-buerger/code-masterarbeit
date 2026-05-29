// Define the base field and number field
Q := RationalField();
R_t<t> := PolynomialRing(Q);
K<zeta> := NumberField(t^4 - t^2 + 1);

i := zeta^3;       // imaginary unit
omega := zeta^4;   // primitive third root of unity

// Define the matrix group G_6
r := Matrix(K, 2, 2, [1, 0, 0, -1]);
r_1 := Matrix(K, 2, 2, [
    (omega/2)*(-1-i), (omega/2)*(1-i),
    (omega/2)*(-1-i), (omega/2)*(-1+i)
]);
G_6 := MatrixGroup<2, K | r, r_1>;

// Define the symplectic reflection group G_6_symp
r_symp := DiagonalJoin(r, Transpose(r^-1));
r_1_symp := DiagonalJoin(r_1, Transpose(r_1^-1));
G_6_symp := MatrixGroup<4, K | r_symp, r_1_symp>;

// Compute the invariants of the symplectic reflection group automatically
IR := InvariantRing(G_6_symp);
invars_auto := FundamentalInvariants(IR);
invars_ideal_1 := ideal<PolynomialRing(IR) | invars_auto>;

// Define the manual invariants as polynomials in the original variables
R<x1, x2, x3, x4> := PolynomialRing(K, 4);
x := [x1, x2, x3, x4];

invar_1 := x1*x3 + x2*x4;
invar_2 := x1^4 + (4*zeta^2 - 2)*x1^2*x2^2 + x2^4;
invar_3 := x3^4 + (-4*zeta^2 + 2)*x3^2*x4^2 + x4^4;
invar_4 := x1^3*x3*x4^2 - x1^2*x2*x4^3 - x1*x2^2*x3^3 + x2^3*x3^2*x4;
invar_5 := x1^6*x4^2 + (4*zeta^2 - 2)*x1^5*x2*x3*x4 - 3*x1^4*x2^2*x3^2 + (-4*zeta^2 + 2)*x1^4*x2^2*x4^2 + 
           4*x1^3*x2^3*x3*x4 + (-4*zeta^2 + 2)*x1^2*x2^4*x3^2 - 3*x1^2*x2^4*x4^2 + (4*zeta^2 - 2)*x1*x2^5*x3*x4 + x2^6*x3^2;
invar_6 := x1^2*x3^6 + 5*x1^2*x3^2*x4^4 + 1/3*(-4*zeta^2 + 2)*x1^2*x4^6 - 2*x1*x2*x3^5*x4 + 
           1/3*(-40*zeta^2 + 20)*x1*x2*x3^3*x4^3 - 2*x1*x2*x3*x4^5 + 1/3*(-4*zeta^2 + 2)*x2^2*x3^6 + 5*x2^2*x3^4*x4^2 + x2^2*x4^6;
invar_7 := x1^9*x3 + (6*zeta^2 - 3)*x1^7*x2^2*x3 + (10*zeta^2 - 5)*x1^6*x2^3*x4 - 9*x1^5*x2^4*x3 - 
           9*x1^4*x2^5*x4 + (10*zeta^2 - 5)*x1^3*x2^6*x3 + (6*zeta^2 - 3)*x1^2*x2^7*x4 + x2^9*x4;
invar_8 := x1*x3^9 - 6*x1*x3^5*x4^4 + (-16*zeta^2 + 8)*x1*x3^3*x4^6 - 3*x1*x3*x4^8 - 3*x2*x3^8*x4 + 
           (-16*zeta^2 + 8)*x2*x3^6*x4^3 - 6*x2*x3^4*x4^5 + x2*x4^9;
invar_9 := x1^12 - 33*x1^8*x2^4 - 33*x1^4*x2^8 + x2^12;
invar_10 := x3^12 - 33*x3^8*x4^4 - 33*x3^4*x4^8 + x4^12;

invars_list := [invar_1, invar_2, invar_3, invar_4, invar_5, invar_6, invar_7, invar_8, invar_9, invar_10];
invars_ideal_2 := ideal<R | invars_list>;

// Compute the relations between the invariants to obtain the presentation using elimination theory.
// Define a big polynomial ring where the first 4 variables will be eliminated
T<tx1,tx2,tx3,tx4, ty1,ty2,ty3,ty4,ty5,ty6,ty7,ty8,ty9,ty10> := PolynomialRing(K, 14, "elim", 4);

// Map the invariants from R into the x-variables of T
map_R_to_T := hom< R -> T | [tx1, tx2, tx3, tx4] >;
invars_in_T := [map_R_to_T(f) : f in invars_list];
ty_vars := [ty1, ty2, ty3, ty4, ty5, ty6, ty7, ty8, ty9, ty10];

// Build the graph ideal: (ty_i - invar_i) using correct bracketed range [1..10]
graph_ideal := ideal< T | [ty_vars[i] - invars_in_T[i] : i in [1..10]] >;

// Eliminate tx1, tx2, tx3, tx4
elim_ideal := EliminationIdeal(graph_ideal, 4);

// Define the target ring S with the desired "grevlexw" ordering working like the ordering available in Oscar
W := [Degree(f) : f in invars_list];
S<y1, y2, y3, y4, y5, y6, y7, y8, y9, y10> := PolynomialRing(K, 10, "grevlexw", W);

// Map the eliminated relations from T into S (sending x-variables to 0)
map_T_to_S := hom< T -> S | [0,0,0,0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10] >;
relations_ideal := ideal< S | [map_T_to_S(f) : f in Basis(elim_ideal)] >;

// Extract the relations basis under the specific ordering
relations_basis := Basis(relations_ideal);
print "Number of relations found:", #relations_basis;

relations_basis_ideal := ideal<S | relations_basis>;