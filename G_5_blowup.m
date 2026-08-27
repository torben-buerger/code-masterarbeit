// Computations needed for the blowup of CC^4/G_5 along the ideal delta_11
// Define the ambient space and some helper functions, these parts are taken from Berry (https://github.com/callumdbberry/partial-resolutions) to be able to reuse some of his methods later on
Prel<Pq> := PolynomialRing(Integers());
Q<q> := NumberField(Pq^2 + 3);
A := AffineSpace(Q, 4);
PQ<x,y,X,Y> := CoordinateRing(A);
q := PQ ! q;

FacNice := function(g);
if g eq 0 then return 0;
elif IsUnit(g) then return 1;
else return Factorisation(g);
end if; end function;

Disp := function(frac);
if frac eq 0 then return 0;
else
    Num := FacNice(Numerator(frac));
    Den := FacNice(Denominator(frac));
    if Type(Den) eq 'RngIntElt' then return Num;
         elif Type(Num) eq 'RngIntElt' then return [[<1,1>], Den];
         else return [Num, Den]; end if;
end if; end function;

rho := hom<PQ -> PQ | x, y, y, -x>;
rhoQ := function(frac);
return rho(Numerator(frac)) / rho(Denominator(frac));
end function;
rhoQF := function(frac);
return Disp(rhoQ(frac));
end function;

phiC := hom<PQ -> PQ | X, Y, x, y>;
phiCQ := function(frac);
return phiC(Numerator(frac)) / phiC(Denominator(frac));
end function;

Bideg := function(A);
if A eq 0 then return [0,0];
else
    Num := PQ ! Numerator(A);
    NumDeg := Exponents(LeadingMonomial(Num));
    Den := PQ ! Denominator(A);
    DenDeg := Exponents(LeadingMonomial(Den));
    return [NumDeg[1] + NumDeg[2] - DenDeg[1] - DenDeg[2],
             NumDeg[3] + NumDeg[4] - DenDeg[3] - DenDeg[4]];
end if; end function;

// Define the polynomials from which the generators of the invariant ring can be computed
F11 := x*X + Y*y;
F40 := x^4 - 2*q*x^2*y^2 + y^4;
F31 := x^3*Y - q*x^2*y*X + q*x*y^2*Y - y^3*X;
F13 := x*Y^3 + q*y*X*Y^2 - q*x*X^2*Y - y*X^3;
F04 := X^4 + 2*q*X^2*Y^2 + Y^4;
F60 := x^5*y - x*y^5;
F33 := x^3*X*Y^2 - x^2*y*Y^3 - x*y^2*X^3 + y^3*X^2*Y;
F06 := X^5*Y - X*Y^5;

F40bar := x^4 + 2*q*x^2*y^2 + y^4;
F31bar := x^3*Y + q*x^2*y*X - q*x*y^2*Y - y^3*X;
F13bar := x*Y^3 - q*y*X*Y^2 + q*x*X^2*Y - y*X^3;
F04bar := X^4 - 2*q*X^2*Y^2 + Y^4;
G22 := x^2*X^2 + q*x^2*Y^2 + q*y^2*X^2 + y^2*Y^2 - 4*x*y*X*Y;
G22bar := x^2*X^2 - q*x^2*Y^2 - q*y^2*X^2 + y^2*Y^2 - 4*x*y*X*Y;
G51 := x^5*X - 5*x^4*y*Y - 5*x*y^4*X + y^5*Y;
G42 := x^4*X*Y - 2*x^3*y*Y^2 + 2*x*y^3*X^2 - y^4*X*Y;
G24 := 2*x^2*X*Y^3 + x*y*X^4 - x*y*Y^4 - 2*y^2*X^3*Y;
G15 := x*X^5 - 5*x*X*Y^4 - 5*y*X^4*Y + y*Y^5;

h := F11;
A1 := 6*F60;
A0 := 3*F33;
A01 := 6*F06;
B1 := -F40*F31 - F40bar*F31bar;
B0 := F31*F13 + F31bar*F13bar;
B01 := -F13*F04 - F13bar*F04bar;
C2 := F40^3 + F40bar^3;
C1 := - F31^3 - F31bar^3;
C0 := F40*F13^2 + F40bar*F13bar^2 + 4*F11^3*F33;
C01 := -F13^3 - F13bar^3;
C02 := F04^3 - F04bar^3;

// Generators of the ideal delta_11 computed in G_5_reflections_2.jl. These will be the ideal along which the blowup is computed
R71 := q*h*A1 + 3*B1;
R44 := 6*h^4 + 4*q*h*A0 + 3*B0;
R120 := -q*A1^2 + 3*C2;
R93 := -2*q*h^2*B1 - q*A1*A0 + 3*C1;
R66 := -12*h^3*A0 + q*A1*A01 - 4*q*A0^2 + 9*C0;

// We list the charts having a strict degree in ChartList, as these charts cannot be covered by other charts and thus have to be included in the necessary charts
GenList := [R71, R44, R120, R93, R66];
ChartList := [R71, R44, R120];

// The method to determine the charts that are not needed to cover the blown up space is due to Berry (https://github.com/callumdbberry/partial-resolutions, e.g. in Blow-up_for_G5.txt), see Remark 3.3.4 in my thesis
n := 2;
RemList := GenList;

GenId := ideal<PQ | GenList>;
GenIdPow := GenId^(n-1);
for Gen in RemList do
	if not Gen in ChartList then
		if Gen^n in ideal<PQ | Exclude(RemList, Gen)> * GenIdPow then
			Exclude(~RemList, Gen);
			"Removed generator of bidegree";
			Bideg(Gen);
		else
			"Included generator of bidegree";
			Bideg(Gen);
		end if;
	end if;
end for; "break";

// List the four resulting charts that cover the blown up space
RemList eq [R71, R44, R120, R66];


// For the charts corresponding to generators that have a strict degree, we can use the following method by Berry (https://github.com/callumdbberry/partial-resolutions, e.g. in Blow-up_for_G5.txt), see Remark 3.3.5 in my thesis
G5Gens := [h, A1, A0, A01, B1, B0, B01, C2, C1, C0, C01, C02];
OrigGens := G5Gens;
NewGens := GenList;

n := 2;

Idealnm1 := GenId^(n-1);

ReqGenList := [];
for Chart in ChartList do
	ChartOGens := [];
	OTestI := Idealnm1 * ideal<PQ | Exclude(NewGens, Chart)>;
	for OGen in OrigGens do
		BVal := OGen*Chart^n in OTestI;
		if not BVal then Append(~ChartOGens, OGen); end if;
	end for;

	ChartNGens := [];
	OtherGens := Exclude(NewGens, Chart);
	for NGen in OtherGens do
		BVal := NGen*Chart^(n-1) in Idealnm1 *
		ideal<PQ | Exclude(OtherGens, NGen)>;
		if not BVal then Append(~ChartNGens, NGen); end if;
	end for;
	
	Data := [[Chart], ChartOGens, ChartNGens];
	ReqGenList cat:= [Data];
	"done for"; Bideg(Chart);
end for;


Chart71 := ReqGenList[1];
Chart44 := ReqGenList[2];
Chart120 := ReqGenList[3];
#Chart71[2];
#Chart71[3];
#Chart44[2];
#Chart44[3];
#Chart120[2];
#Chart120[3];

Chart71[2] eq [h, A1, A0, B1, C2];
Chart71[3] eq [R44, R120, R93, R66];

Chart44[2] eq [h, A01, B01, C02];
Chart44[3] eq [R71, R120, R93, R66];

Chart120[2] eq [h, A1, B1, C2];
Chart120[3] eq [R71, R44, R93];

// The chart X66 has a strict degree 6x+y after excluding R44/R66. Therefore we apply the above method to this smaller GenList
OrigGens := G5Gens;
NewGens := Exclude(GenList, R44);
ChartList := [R66];

n := 4;
GenList := Exclude(GenList, R44);
GenId := ideal<PQ | GenList>;

Idealnm1 := GenId^(n-1);

ReqGenList := [];
for Chart in ChartList do
	ChartOGens := [];
	OTestI := Idealnm1 * ideal<PQ | Exclude(NewGens, Chart)>;
	for OGen in OrigGens do
		BVal := OGen*Chart^n in OTestI;
		if not BVal then Append(~ChartOGens, OGen); end if;
	end for;
	
	ChartNGens := [];
	OtherGens := Exclude(NewGens, Chart);
	for NGen in OtherGens do
		BVal := NGen*Chart^(n-1) in Idealnm1 *
		ideal<PQ | Exclude(OtherGens, NGen)>;
		if not BVal then Append(~ChartNGens, NGen); end if;
	end for;
	
	Data := [[Chart], ChartOGens, ChartNGens];
	ReqGenList cat:= [Data];
	"done for"; Bideg(Chart);
end for;

// Resulting generators after taking out R44/R66 to get a strict degree. So for the full generating set, we have to add this generator back
Chart66 := ReqGenList[1];
Chart66[2] eq [h, A0, A01, B01, C01, C02];
Chart66[3] eq [R71, R93];


// From now on we work in the field of fractions associated to PQ because considering the charts leads to rational functions
F_PQ := FieldOfFractions(PQ);

// Transfer the previously computed generators for the coordinate algebras to the set field of fractions
X120_gens := [ F_PQ!h, F_PQ!A1, F_PQ!B1, F_PQ!C2,
                F_PQ!(R71/R120), F_PQ!(R44/R120), F_PQ!(R93/R120) ];

X71_gens := [ F_PQ!h, F_PQ!A1, F_PQ!A0, F_PQ!B1, F_PQ!C2,
              F_PQ!(R44/R71), F_PQ!(R120/R71), F_PQ!(R93/R71), F_PQ!(R66/R71) ];

X44_gens := [ F_PQ!h, F_PQ!A01, F_PQ!B01, F_PQ!C02,
               F_PQ!(R71/R44), F_PQ!(R120/R44), F_PQ!(R93/R44), F_PQ!(R66/R44) ];

X66_gens := [ F_PQ!h, F_PQ!A0, F_PQ!A01, F_PQ!B01, F_PQ!C01, F_PQ!C02,
               F_PQ!(R71/R66), F_PQ!(R44/R66), F_PQ!(R93/R66) ];

FindChartRelations := function(gens, chart_name : MaxDegree := 0)
    k := #gens;
    F_PQ := FieldOfFractions(PQ);
    gens_frac := [F_PQ ! g : g in gens];
    
    // Create the relation ring: Q[Z_1, ..., Z_k] with variable names suitable to the charts
    RelRing := PolynomialRing(Q, k);
    AssignNames(~RelRing, [ chart_name * "_" * IntegerToString(i) : i in [1..k] ]);
    
    // First, we determine all monomials in the relation ring up to the given degree bound
    mons := &cat[ Setseq(MonomialsOfDegree(RelRing, d)) : d in [0..MaxDegree] ];
    
    // We evaluate each monomial at the generators of the coordinate algebra for which we want to find the relations
    evals := [];
    for m in mons do
        exp := Exponents(m);
        val := F_PQ ! 1;
        for i in [1..k] do
            if exp[i] gt 0 then 
                val *:= gens_frac[i]^exp[i]; 
            end if;
        end for;
        Append(~evals, val);
    end for;
    
    // Clear denominators across all evaluated rational functions
    LCM_den := PQ ! 1;
    for v in evals do
        if v ne 0 then
            LCM_den := Lcm(LCM_den, Denominator(v));
        end if;
    end for;
    
    // Using the polynomials with cleared denominators, we work in the field original ring to collect all x, y, X, Y monomials
    polys := [ PQ ! (v * LCM_den) : v in evals ];
    all_x_mons_set := {};
    for p in polys do
        for m in Monomials(p) do
            Include(~all_x_mons_set, m);
        end for;
    end for;
    all_x_mons := Setseq(all_x_mons_set);
    
    if #all_x_mons eq 0 then  // in this case, no relations can be found
        return ideal<RelRing | >; 
    end if;
    
    // The coefficient matrix is index as follows: rows  are given by the monomials in RelRing, and columns by the above monomials in PQ
    // Its kernel yields to coefficients of monimials such that the associated linear combination of monomial becomes 0 and thus forms a relation
    mat := Matrix(Q, [[ MonomialCoefficient(p, m) : m in all_x_mons ] : p in polys]);
    K := Nullspace(mat);
    
    relations := [];
    for vec in Basis(K) do
        rel := RelRing ! 0;
        for i in [1..#mons] do
            rel +:= vec[i] * mons[i];
        end for;
        Append(~relations, rel);
    end for;
    
    I_rel := ideal<RelRing | relations>;
    return ideal<RelRing | Basis(I_rel)>;
end function;

Ideal_X120 := FindChartRelations(X120_gens, "X120" : MaxDegree := 5);
Dimension(Ideal_X120) eq 4;
IsPrime(Ideal_X120);

Ideal_X71 := FindChartRelations(X71_gens,  "X71"  : MaxDegree := 5);
Dimension(Ideal_X71) eq 4;
IsPrime(Ideal_X71);

Ideal_X44 := FindChartRelations(X44_gens,  "X44"  : MaxDegree := 7);
Dimension(Ideal_X44) eq 4;

Ideal_X66 := FindChartRelations(X66_gens,  "X66"  : MaxDegree := 6);
Dimension(Ideal_X66) eq 4;

// As the primality tests for X44 and X66 are too expensive in characteristic 0, we pass to a large prime field containing a root of -3 to check primality
IsPrimeModular := function(I_rel)
    // This function can be easily adapted for number fields Q(sqrt(z)) for an integer z. 
    // First, we have to find a suitable prime field and define the corresponding modular relation ring and pass the relations to it
    p := 10000;
    repeat
        p := NextPrime(p);
    until (p mod 3 eq 1) and IsSquare(GF(p) ! -3);
    
    Fp := GF(p);
    q_val := SquareRoot(Fp ! -3);
    
    RelRing := Generic(I_rel);
    k := Rank(RelRing);
    RelRing_p := PolynomialRing(Fp, k);
    AssignNames(~RelRing_p, [ "Z_" * IntegerToString(i) : i in [1..k] ]);
    
    // Helper function for the transition of q3 = sqrt(-3) to its counterpart in the prime field
    Mapq3ToFp := function(c)
        coeffs := ElementToSequence(c); // Decompose a coefficient of Q(sqrt(-3)) of the form c = a + b*q3 and transfer both parts individually to Fp
        a := coeffs[1];
        a_fp := Fp ! Numerator(a) / Fp ! Denominator(a);
        if #coeffs eq 1 then
            return a_fp;
        else
            b := coeffs[2];
            b_fp := Fp ! Numerator(b) / Fp ! Denominator(b);
            return a_fp + b_fp * q_val;
        end if;
    end function;
    
    // We compute the reduction of the generators of the characteristic 0 ideal to the modular relation ring
    gens_p := [];
    for f in Basis(I_rel) do
        f_p := RelRing_p ! 0;
        for m in Monomials(f) do
            coeff := MonomialCoefficient(f, m);
            exp := Exponents(m);
            mon_p := &*[ RelRing_p.i^exp[i] : i in [1..k] ];
            f_p +:= Mapq3ToFp(coeff) * mon_p;
        end for;
        Append(~gens_p, f_p);
    end for;
    
    // Define the ideal and check its primality
    I_p := ideal< RelRing_p | gens_p >;
    res := IsPrime(I_p);
    
    if res then
        print "[=== SUCCESS ===] Ideal is prime over GF(p)! Relation set is complete.";
    else
        print "[-] FAIL: Ideal is not prime over GF(p). Higher degree relations needed.";
    end if;
    return res;
end function;
IsPrimeModular(Ideal_X44);
IsPrimeModular(Ideal_X66);

// Redefine the relation rings with an alternative polynomial ordering to possibly facilitate the generating sets
X71_revlex<x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9> := PolynomialRing(Q, 9, "grevlex");
transfer_X71 := hom<Parent(Basis(Ideal_X71)[1]) -> X71_revlex | x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9 >;
Ideal_X71_new := ideal< X71_revlex | [transfer_X71(Basis(Ideal_X71)[i]) : i in [1..#Basis(Ideal_X71)]] >;

X120_revlex<x_1, x_2, x_3, x_4, x_5, x_6, x_7> := PolynomialRing(Q, 7, "grevlex");
transfer_X120 := hom<Parent(Basis(Ideal_X120)[1]) -> X120_revlex | x_1, x_2, x_3, x_4, x_5, x_6, x_7 >;
Ideal_X120_new := ideal< X120_revlex | [transfer_X120(Basis(Ideal_X120)[i]) : i in [1..#Basis(Ideal_X120)]] >;

X44_revlex<x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8> := PolynomialRing(Q, 8, "grevlex");
transfer_X44 := hom<Parent(Basis(Ideal_X44)[1]) -> X44_revlex | x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8 >;
Ideal_X44_new := ideal< X44_revlex | [transfer_X44(Basis(Ideal_X44)[i]) : i in [1..#Basis(Ideal_X44)]] >;

X66_revlex<x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9> := PolynomialRing(Q, 9, "grevlex");
transfer_X66 := hom<Parent(Basis(Ideal_X66)[1]) -> X66_revlex | x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9 >;
Ideal_X66_new := ideal< X66_revlex | [transfer_X66(Basis(Ideal_X66)[i]) : i in [1..#Basis(Ideal_X66)]] >;

// Heurstic function usable to find smaller generating sets for an ideal by trying out if
function ReduceGenerators(I)
    G := GroebnerBasis(I);
    G := Sort(G, func<a,b | TotalDegree(b) - TotalDegree(a)>); // start with the generators of the largest degree
    keep := G;
    for g in G do
        rest := Exclude(keep, g);
        if #rest eq 0 then continue; end if;
        J := ideal<Universe(rest) | rest>;
        if g in J then          // test if g is already contained in the ideal generated by the remaining generators, if yes it is superfluous and removed
            keep := rest;
        end if;
    end for;
    return keep;
end function;

// As for the OSCAR computations (e.g. blowupsW1.jl, blowupsW2.jl) we reduce the embedding dimension. This function is a direct translation of the OSCAR versoin
MinimalEmbedding := function(I)
    R := Generic(I);
    nvars := Rank(R);
    ring_vars := [R.i : i in [1..nvars]];
    active_vars := ring_vars;
    current_gens := Basis(I);
    
    // Eliminated variables and their substitution expressions are store as tuples of the form <Var, Expr>
    eliminated_pairs := [];
    
    cleaning := true;
    while cleaning do
        cleaning := false;
        
        for g in current_gens do
            if g eq 0 then continue; end if;
            
            found_v := 0;
            found_idx := 0;
            subst_expr := R ! 0;
            
            for v in active_vars do
                if Degree(g, v) eq 1 then
                    c := Derivative(g, v);
                    
                    // Check if c is a non-zero constant
                    if TotalDegree(c) eq 0 and c ne 0 then
                        rest := g - c*v;
                        
                        // Extract constant scalar to divide it in the base field Q
                        c_val := MonomialCoefficient(c, R ! 1);
                        subst_expr := -rest * (c_val^-1);
                        
                        found_v := v;
                        found_idx := Index(ring_vars, v);
                        break;
                    end if;
                end if;
            end for;
            
            if found_idx gt 0 then
                printf "Eliminating variable: %o\n", ring_vars[found_idx];
                
                Append(~eliminated_pairs, <ring_vars[found_idx], subst_expr>);
                Exclude(~active_vars, found_v);
                
                // Substitute the expression for the variable to be eliminated across the remaining generators
                new_gens := [];
                for poly in current_gens do
                    p_sub := Evaluate(poly, found_idx, subst_expr);
                    if p_sub ne 0 then
                        Append(~new_gens, p_sub);
                    end if;
                end for;
                current_gens := new_gens;
                
                cleaning := true;
                break;
            end if;
        end for;
    end while;
    
    // As in OSCAR, we return the remaining variables, the ideal after the substitutions and the elimination rules
    return active_vars, ideal<R | current_gens>, eliminated_pairs;
end function;

// For all charts, we do the same computations: 1. find a smaller generating set, 2. determine a minimal embedding and a smaller generating therof
// 3. define the corresponding affine scheme and compute its singular locus and find a smaller generating set thereof, 4. if applicable verify the branching results in Section 3.2 by elimination

// Chart X_71
Red_X71 := ReduceGenerators(Ideal_X71_new);
Ideal_Red_X71 := ideal< X71_revlex | Red_X71 >;
Ideal_Red_X71 eq Ideal_X71_new;  // returns true
vars_X71, min_ideal_X71, elim_X71 := MinimalEmbedding(Ideal_Red_X71);
min_Red_X71 := ReduceGenerators(min_ideal_X71); 
min_Red_ideal_X71 := ideal< X71_revlex | min_Red_X71 >;
min_Red_ideal_X71 eq min_ideal_X71;  // returns true
// Affine scheme and singular locus
A_X71 := AffineSpace(X71_revlex);
X71 := Scheme(A_X71, min_Red_ideal_X71);
Sing_X71 := SingularSubscheme(X71);
I_sing_X71 := Ideal(Sing_X71);
I_sing_X71_reduced := Radical(I_sing_X71);
Red_I_sing_X71_reduced := ReduceGenerators(I_sing_X71_reduced);
#Red_I_sing_X71_reduced;
I_red_sing_X71 := ideal<X71_revlex|Red_I_sing_X71_reduced>;
// Branching and elimination
Branch1_ideal := ideal<X71_revlex|X71_revlex.7, X71_revlex.9, X71_revlex.8 + X71_revlex!(1/q)*(X71_revlex.1)^2, X71_revlex.4 - X71_revlex!(1/q)*(X71_revlex.1)*(X71_revlex.2)>;
Branch1 := EliminationIdeal(I_red_sing_X71 + Branch1_ideal, {X71_revlex.1, X71_revlex.2, X71_revlex.6});
Branch2_ideal := ideal<X71_revlex|X71_revlex.2 + X71_revlex.1*X71_revlex.7, X71_revlex.8 + X71_revlex.6*X71_revlex.7/4>;
Branch2 := EliminationIdeal(I_red_sing_X71 + Branch2_ideal, {X71_revlex.1, X71_revlex.6, X71_revlex.7});


// Chart X_120
Red_X120 := ReduceGenerators(Ideal_X120_new);
Ideal_Red_X120 := ideal< X120_revlex | Red_X120 >;
Ideal_Red_X120 eq Ideal_X120_new;   // returns true 
vars_X120, min_ideal_X120, elim_X120 := MinimalEmbedding(Ideal_Red_X120); 
min_Red_X120 := ReduceGenerators(min_ideal_X120); 
min_Red_ideal_X120 := ideal< X120_revlex | min_Red_X120 >;
min_Red_ideal_X120 eq min_ideal_X120;  // returns true
// Affine scheme and singular locus
A_X120 := AffineSpace(X120_revlex);
X120 := Scheme(A_X120, min_Red_ideal_X120);
Sing_X120 := SingularSubscheme(X120);
I_sing_X120 := Ideal(Sing_X120);
I_sing_X120_reduced := Radical(I_sing_X120);
Red_I_sing_X120_reduced := ReduceGenerators(I_sing_X120_reduced);
#Red_I_sing_X120_reduced;
I_red_sing_X120 := ideal<X120_revlex|Red_I_sing_X120_reduced>;
Dimension(I_red_sing_X120);
IsPrime(I_red_sing_X120);
// Elimination
Elim_X120 := EliminationIdeal(I_red_sing_X120, {X120_revlex.2, X120_revlex.5});


// Chart X_44
Red_X44 := ReduceGenerators(Ideal_X44_new);
Ideal_Red_X44 := ideal< X44_revlex | Red_X44 >;
Ideal_Red_X44 eq Ideal_X44_new; // returns true 
vars_X44, min_ideal_X44, elim_X44 := MinimalEmbedding(Ideal_Red_X44); 
min_Red_X44 := ReduceGenerators(min_ideal_X44);
min_Red_ideal_X44 := ideal< X44_revlex | min_Red_X44 >;
min_Red_ideal_X44 eq min_ideal_X44;  // returns true
// Affine scheme and singular locus
A_X44 := AffineSpace(X44_revlex);
X44 := Scheme(A_X44, min_Red_ideal_X44);
Sing_X44 := SingularSubscheme(X44);
I_sing_X44 := Ideal(Sing_X44);
I_sing_X44_reduced := Radical(I_sing_X44);
Red_I_sing_X44_reduced := ReduceGenerators(I_sing_X44_reduced);
#Red_I_sing_X44_reduced;
I_red_sing_X44 := ideal<X44_revlex|Red_I_sing_X44_reduced>;
// Branching and elimination
Branch11_ideal := ideal<X44_revlex|X44_revlex.8, X44_revlex.3 - X44_revlex.1*X44_revlex.2*X44_revlex!(1/q)>;
Brach11 := EliminationIdeal(I_red_sing_X44 + Branch11_ideal, {X44_revlex.1, X44_revlex.2,X44_revlex.5});
Branch2_ideal := ideal<X44_revlex|X44_revlex.8 + 3*(X44_revlex.1)^2*X44_revlex!(1/q), X44_revlex.3 + X44_revlex.1*X44_revlex.2*X44_revlex!(1/q)>;
Branch2 := EliminationIdeal(I_red_sing_X44 + Branch2_ideal, {X44_revlex.1, X44_revlex.2,X44_revlex.5});


// Chart X_66
Red_X66 := ReduceGenerators(Ideal_X66_new); // very slow for degree 6 relations, computing the minimal embedding first works better
Ideal_Red_X66 := ideal< X66_revlex | Red_X66 >;
Ideal_Red_X66 eq Ideal_X66_new;
vars_X66, min_ideal_X66, elim_X66 := MinimalEmbedding(Ideal_Red_X66); 
min_Red_X66 := ReduceGenerators(min_ideal_X66);
min_Red_ideal_X66 := ideal< X66_revlex | min_Red_X66 >;
min_Red_ideal_X66 eq min_ideal_X66;  // returns true
// Affine scheme and singular locus
A_X66 := AffineSpace(X66_revlex);
X66 := Scheme(A_X66, min_Red_ideal_X66);
Sing_X66 := SingularSubscheme(X66);
I_sing_X66 := Ideal(Sing_X66);
I_sing_X66_reduced := Radical(I_sing_X66);
Red_I_sing_X66_reduced := ReduceGenerators(I_sing_X66_reduced);
#Red_I_sing_X66_reduced;

