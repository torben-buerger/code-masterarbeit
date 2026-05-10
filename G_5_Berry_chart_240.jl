# This file only contains the computations needed to establish the relations of the variables needed to show that the chart X_24,0 is smoot
using Oscar
Q_t, t = polynomial_ring(QQ, :t);
K, zeta = number_field(t^8 - t^4 + 1, "zeta");
i = zeta^6;      # imaginary unit
sigma = zeta^8;  # primitive third root of unity
q2 = zeta^3 + zeta^(-3); # sqrt(2)
q3 = zeta^8 - zeta^(-8); # sqrt(3)

# We define the previously computed semi-invariants (eigenvectors of c and c_conj), the resulting invariants and the generators of the ideal in which the blowup is comoputed
R, x = polynomial_ring(K, :x => (1:4));

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

# Define a function that mimics the complex conjugation
f_conj = hom(K, K, -zeta^7+zeta^3);
function conj_poly(poly)
    return map_coefficients(f_conj, poly)
end

F40_conj = conj_poly(F40);
F31_conj = conj_poly(F31);
G22_conj = conj_poly(G22);
F13_conj = conj_poly(F13);
F04_conj = conj_poly(F04);

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

# The blowup will be computed in the ideal corresponding to the singular locus which is given by the following semi invariants
R44 = G22*G22_conj;
R93 = F31*F40_conj*G22_conj;
R93prime = F40*F31_conj*G22;
R66 = F04*F40_conj*G22_conj;
R66prime = F40*F04_conj*G22;
R39 = F04*F13_conj*G22_conj;
R39prime = F13*F04_conj*G22;
R131 = F40*F40_conj*G51;
R113 = F04*F04_conj*G15;
R142 = F40*F31*F40_conj*F31_conj;
R214 = F13*F04*F13_conj*F04_conj;
R191 = F40*F31*F40_conj^3;
R119 = F13*F04*F04_conj^3;
R240 = F40^3*F40_conj^3;
R024 = F04^3*F04_conj^3;

# As we aim to find relations after assuming R_24,0 != 0, by expressing all variables using A1, C2, R131/R240, and R191/R240, we have to find expressions for all other variables in these four generators
# To do so we create a list of all monomials in the four generators that have the matching bidegree of the investigated polynomial. Note that we do not work in the localized ring, but give expressions that have to be divided by a suitable power of R240
# This functions enumerates all the monomials of a given bidgree

function get_bidegree_monomials(deg_1::Int, deg_2::Int)    
    monomials = elem_type(R)[]
    # Iterate through possible R191 exponents (d)
    for d in 0:deg_2
        c = deg_2 - d
        r = deg_1 + 11*c + 5*d
        if mod(r, 6) != 0  # We have to solve 6*a+12*b = r, thus a solution only exists if r is divisible by 6
            return []
        end
        r0 = div(r, 6)
        # If r0 is negative, no non-negative a, b exist for this c
        if r0 >= 0
            # Iterate through possible C2 exponents (b)
            for b in 0:div(r0, 2)
                a = r0 - 2*b
                d = deg_2 - c
                push!(monomials, A1^a*C2^b*R131^c*R191^d)
            end
        end
    end
    return monomials
end

# 1. h of bidegree (1,1)
monomials_h = get_bidegree_monomials(1,1);
ideal_h = ideal(R, monomials_h);
print(ideal_membership(h*R240, ideal_h));
coeffs_h = coordinates(h*R240, ideal_h);
print((monomials_h*transpose(coeffs_h))[1] == h*R240);

# 2. A0 of bidegree (3,3)
monomials_A0 = get_bidegree_monomials(3,3);
ideal_A0 = ideal(R, monomials_A0);
print(ideal_membership(A0*R240^3, ideal_A0));
coeffs_A0 = coordinates(A0*R240^3, ideal_A0);
print((monomials_A0*transpose(coeffs_A0))[1] == A0*R240^3);

# 3. Aneg1 of bidegree (0,6)
monomials_Aneg1 = get_bidegree_monomials(0,6);
ideal_Aneg1 = ideal(R, monomials_Aneg1);
print(ideal_membership(Aneg1*R240^6, ideal_Aneg1));
coeffs_Aneg1 = coordinates(Aneg1*R240^6, ideal_Aneg1);
print((monomials_Aneg1*transpose(coeffs_Aneg1))[1] == Aneg1*R240^6);

# 4. B1 of bidegree (7,1)
monomials_B1 = get_bidegree_monomials(7,1);
ideal_B1 = ideal(R, monomials_B1);
print(ideal_membership(B1*R240, ideal_B1));
coeffs_B1 = coordinates(B1*R240, ideal_B1);
print((monomials_B1*transpose(coeffs_B1))[1] == B1*R240);

# 5. B0 of bidegree (4,4)
monomials_B0 = get_bidegree_monomials(4,4);
ideal_B0 = ideal(R, monomials_B0);
print(ideal_membership(B0*R240^4, ideal_B0));
coeffs_B0 = coordinates(B0*R240^4, ideal_B0);
print((monomials_B0*transpose(coeffs_B0))[1] == B0*R240^4);

# 6. Bneg1 of bidegree (1,7)
monomials_Bneg1 = get_bidegree_monomials(1,7);
ideal_Bneg1 = ideal(R, monomials_Bneg1);
print(ideal_membership(Bneg1*R240^7, ideal_Bneg1));
coeffs_Bneg1 = coordinates(Bneg1*R240^7, ideal_Bneg1);
print((monomials_Bneg1*transpose(coeffs_Bneg1))[1] == Bneg1*R240^7);

# 7. C1 of bidegree (9,3)
monomials_C1 = get_bidegree_monomials(9,3);
ideal_C1 = ideal(R, monomials_C1);
print(ideal_membership(C1*R240^3, ideal_C1));
coeffs_C1 = coordinates(C1*R240^3, ideal_C1);
print((monomials_C1*transpose(coeffs_C1))[1] == C1*R240^3);

# 8. C0 of bidegree (6,6)
monomials_C0 = get_bidegree_monomials(6,6);
ideal_C0 = ideal(R, monomials_C0);
print(ideal_membership(C0*R240^6, ideal_C0));
coeffs_C0 = coordinates(C0*R240^6, ideal_C0);
print((monomials_C0*transpose(coeffs_C0))[1] == C0*R240^6);

# 9. Cneg1 of bidegree (3,9)
monomials_Cneg1 = get_bidegree_monomials(3,9);
ideal_Cneg1 = ideal(R, monomials_Cneg1);
print(ideal_membership(Cneg1*R240^9, ideal_Cneg1));
coeffs_Cneg1 = coordinates(Cneg1*R240^9, ideal_Cneg1);
print((monomials_Cneg1*transpose(coeffs_Cneg1))[1] == Cneg1*R240^9);

# 10. Cneg2 of bidegree (0,12)
monomials_Cneg2 = get_bidegree_monomials(0,12);
ideal_Cneg2 = ideal(R, monomials_Cneg2);
print(ideal_membership(Cneg2*R240^12, ideal_Cneg2));
coeffs_Cneg2 = coordinates(Cneg2*R240^12, ideal_Cneg2);
print((monomials_Cneg2*transpose(coeffs_Cneg2))[1] == Cneg2*R240^12);

# 11. R44/R240 of bidegree (-20,4)
monomials_R44 = get_bidegree_monomials(-20,4);
ideal_R44 = ideal(R, monomials_R44);
print(ideal_membership(R44*R240^3, ideal_R44));
coeffs_R44 = coordinates(R44*R240^3, ideal_R44);
print((monomials_R44*transpose(coeffs_R44))[1] == R44*R240^3);

# 12. R93/R240 of bidegree (-15,3)
monomials_R93 = get_bidegree_monomials(-15,3);
ideal_R93 = ideal(R, monomials_R93);
print(ideal_membership(R93*R240^2, ideal_R93));
coeffs_R93 = coordinates(R93*R240^2, ideal_R93);
print((monomials_R93*transpose(coeffs_R93))[1] == R93*R240^2);

# 13. R93prime/R240 of bidegree (-15,3)
monomials_R93prime = get_bidegree_monomials(-15,3);
ideal_R93prime = ideal(R, monomials_R93prime);
print(ideal_membership(R93prime*R240^2, ideal_R93prime));
coeffs_R93prime = coordinates(R93prime*R240^2, ideal_R93prime);
print((monomials_R93prime*transpose(coeffs_R93prime))[1] == R93prime*R240^2);

# 14. R66/R240 of bidegree (-18,6)
monomials_R66 = get_bidegree_monomials(-18,6);
ideal_R66 = ideal(R, monomials_R66);
print(ideal_membership(R66*R240^5, ideal_R66));
coeffs_R66 = coordinates(R66*R240^5, ideal_R66);
print((monomials_R66*transpose(coeffs_R66))[1] == R66*R240^5);

# 15. R66prime/R240 of bidegree (-18,6)
monomials_R66prime = get_bidegree_monomials(-18,6);
ideal_R66prime = ideal(R, monomials_R66prime);
print(ideal_membership(R66prime*R240^5, ideal_R66prime));
coeffs_R66prime = coordinates(R66prime*R240^5, ideal_R66prime);
print((monomials_R66prime*transpose(coeffs_R66prime))[1] == R66prime*R240^5);

# 16. R39/R240 of bidegree (-21,9)
monomials_R39 = get_bidegree_monomials(-21,9);
ideal_R39 = ideal(R, monomials_R39);
print(ideal_membership(R39*R240^8, ideal_R39));
coeffs_R39 = coordinates(R39*R240^8, ideal_R39);
print((monomials_R39*transpose(coeffs_R39))[1] == R39*R240^8);

# 17. R39prime/R240 of bidegree (-21,9)
monomials_R39prime = get_bidegree_monomials(-21,9);
ideal_R39prime = ideal(R, monomials_R39prime);
print(ideal_membership(R39prime*R240^8, ideal_R39prime));
coeffs_R39prime = coordinates(R39prime*R240^8, ideal_R39prime);
print((monomials_R39prime*transpose(coeffs_R39prime))[1] == R39prime*R240^8);

# 18. R113/R240 of bidegree (-23,13)
monomials_R113 = get_bidegree_monomials(-23,13);
ideal_R113 = ideal(R, monomials_R113);
print(ideal_membership(R113*R240^12, ideal_R113));
coeffs_R113 = coordinates(R113*R240^12, ideal_R113);
print((monomials_R113*transpose(coeffs_R113))[1] == R113*R240^12);

# 19. R142/R240 of bidegree (-10,2)
monomials_R142 = get_bidegree_monomials(-10,2);
ideal_R142 = ideal(R, monomials_R142);
print(ideal_membership(R142*R240^1, ideal_R142));
coeffs_R142 = coordinates(R142*R240^1, ideal_R142);
print((monomials_R142*transpose(coeffs_R142))[1] == R142*R240^1);

# 20. R214/R240 of bidegree (-22,14)
monomials_R214 = get_bidegree_monomials(-22,14);
ideal_R214 = ideal(R, monomials_R214);
print(ideal_membership(R214*R240^13, ideal_R214));
coeffs_R214 = coordinates(R214*R240^13, ideal_R214);
print((monomials_R214*transpose(coeffs_R214))[1] == R214*R240^13);

# 21. R119/R240 of bidegree (-23,19)
monomials_R119 = get_bidegree_monomials(-23,19);
ideal_R119 = ideal(R, monomials_R119);
print(ideal_membership(R119*R240^18, ideal_R119));
coeffs_R119 = coordinates(R119*R240^18, ideal_R119);
print((monomials_R119*transpose(coeffs_R119))[1] == R119*R240^18);

# 22. R024/R240 of bidegree (-24,24)
monomials_R024 = get_bidegree_monomials(-24,24);
ideal_R024 = ideal(R, monomials_R024);
print(ideal_membership(R024*R240^23, ideal_R024));
coeffs_R024 = coordinates(R024*R240^23, ideal_R024);
print((monomials_R024*transpose(coeffs_R024))[1] == R024*R240^23);
