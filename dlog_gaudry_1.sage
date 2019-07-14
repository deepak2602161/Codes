import time
from sage.rings.polynomial.polydict import PolyDict

def dict_to_poly(rep, R):
    x = R.gens()
    f = 0
    for t in rep:
        g = 1
        for j in range(len(t)-1):
            g = g*(x[j]**t[j])
        g = g*t[-1]
        f += g
    return f
def curve(i, n, A,  B):
    p = next_prime(i)
    F = GF(p^4, 'a')
    E = EllipticCurve(F, [A, B])
    return [p, E, F]

# this function works for functions defined over multivariate polynomial ring in x, y over a given field.
def weil_descent(f1, F):
    L = F.vector_space()
    R2.<x, y> = PolynomialRing(F, 2, order = 'lex')
    f2 = f1(x, y)
    R3.<y> = PolynomialRing(F)
    R4.<x> = PolynomialRing(R3)
    f = f2
    print(f.parent())
    cf_x = f.coefficients()
    cf_y = []
    for i in range(len(cf_x)):
        cf_y = cf_y + [cf_x[i].coefficients()]
    alpha = []
    beta = []
    
    for i in range(len(cf_y)):
       for j in range(len(cf_y[i])):
           c = L.coordinates(cf_y[i][j])
           c = [k*x^i*y^j for k in c]
           alpha = alpha + c

    k = 0
    n = F.degree()
    for j in range(n):
        for i in range(len(alpha)/n):
            k = k + alpha[n*i + j]
        beta = beta + [k] 
        k = 0 
    return beta

# this function works for any general function defined over multivariate polynomial ring in n variables over a given field.
def weil_descent_1(f, m, F):
    R = PolynomialRing(F, 'x', m, order = 'lex')
    L = F.vector_space()
    f1 = f(R.gens())
    mon_list = f1.monomials()
    coef_list = f1.coefficients()
    alpha = []
    for i in range(len(mon_list)):
        c0 = mon_list[i]
        c1 = coef_list[i]
        c2 = L(c1)
        c2 = [c0*j for j in c2]
        alpha += c2

    beta = []
    k = 0
    n = F.degree()
    for j in range(n):
        for i in range(len(alpha)/n):
            k = k + alpha[n*i + j]
        beta = beta + [k] 
        k = 0 
    return beta


def summation_poly(k, R, A, B):
    F = R.base_ring()
    a = F.gens()[0]
    x = R.gens()
    S = [1, x[0] - x[1]]
    S = S + [(x[0] - x[1])^2*x[2]^2 - 2*((x[0] + x[1])*(x[0]*x[1] + A) + 2*B)*x[2] + (x[0]*x[1] - A)^2 - 4*B*(x[0] + x[1])]
    S1 = PolynomialRing(R, 'X')
    X = S1.gens()[0]
    t1 = list(x) 
    t2 = list(x)
    t1[k-2] = X 
    t2 = [x[k-2], x[k-1], X] + t2[3:]
    if k==3:
        return S[k-1]
    elif k == 4:
        f_1 = S[2](t1)
        f_2 = S[2](t2)
        M1 = f_1.sylvester_matrix(f_2, X)
        S = S + [M1.determinant()]
        return S[k-1]
    else:
        f = summation_poly(k-1, R, A, B)
        S = S + [f]
        M = f(t1).sylvester_matrix(S[2](t2), X) 
        S = S + [M.determinant()]
    return S[-1]

# G = grobner basis of I, I is an ideal over R, we are computing ith elimination ideal(I_i), 
# where I_i = G intersection F[x_i, ..., x_(n-1)]), R = F[x_0, x_1, ..., x_(m-1)], F = GF(p^n).
def elimination_ideal(G, R, i):
    G1 = []
    for j in range(len(G)):
        count = 0
        for k in range(i):
            if G[j].degree(R.gens()[k]) >= 1:
                break
            else:
                count += 1
        if count == i:
            G1 = G1 + [G[j]]
    return G1   

def solve_groebner_basis(S, R):
    solution = []
    n = len(R.gens())
    x = R.gens()
    t = list(x)
    F = R.base_ring()
    a = F.gens()[0]
    R1.<X> = F[]
    I = R.ideal(S)
    basis = I.groebner_basis()
    print(basis)
    if 1 in basis:
        return solution
    else:    
        J = elimination_ideal(basis, R, n-1)
        print(J)
        if J == []:
            s = (F.characteristic())**(F.degree())
            if F.degree() == 1:
                solution = [[i*a] for i in range(s)]
            else:
                solution = [[0]] + [[a**i] for i in range(s-1)]        
            print(solution)
        else:
            t[:(n-1)] = [0]*(n-1)
            t[-1] = X
            f = J[0](t)    
            roots = [[y[0]] for y in f.roots()]
            solution += roots
            print(solution)
        for i in range(n-1)[::-1]:
            if solution == []:
                return solution
                break
            K = elimination_ideal(basis, R, i+1)
            J = elimination_ideal(basis, R, i)
            # print(J)
            solution1 = []
            if J == K:
                for j in solution:
                    solution1 = [[0] + j] + [[a**i] + j for i in range(s-1)]
            else:        
                for j in solution:
                    l = []
                    for k in J:
                        t1 = list(x)
                        t1[:i] = [0]*i
                        t1[i] = X 
                        t1[i+1:] = j 
                        l += [k(t1)]
                    # print(l)
                    I1 = R1.ideal(l)
                    g = I1.groebner_basis()[0]
                    # print(g)
                    if g == 1:
                        break        
                    elif g.roots() == []:
                        break
                    else:    
                        roots = [[y[0]] for y in g.roots()]
                        for s in roots:
                            solution1 += [s + j]
                            # print(solution1)
            solution = solution1
            # print(solution)
    return solution      

def factor_base(F, A, B):
    E = EllipticCurve(F, [A, B])
    n = F.degree()
    p = F.characteristic()
    a = F.gens()[0]
    F_base = F.base_ring()
    factor_base = []
    for i in range(p):
        t = F(i**3 + A*x + B)
        if F(t).nth_root(2) 




# eliptic curve discrete log computation using index calculus algorithm given by Gaudry.
def relations(P, Q, F, k, A, B):
    count = 0 #number of relations
    relations = []
    E = EllipticCurve(F, [A, B])
    n = F.degree()
    p = F.characteristic()
    a = F.gens()[0]
    R = PolynomialRing(F, 'x', k, order = 'lex')
    f = summation_poly(k, R, A, B)
    card = E.cardinality()
    while (count <= p):



#F = F_{q^n2}, q = p^n1(Finite Field) or number field or complex field.    
n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
N = input("insert the size of the prime.")
A, B = input("insert the parameters of required elliptic curve.")
n = input("number of variables in the polynomial ring.")
[p, E, F] = curve(N, n, A, B)
Fn = GF(p^n2, 'a')
F = Fn.extension(n1, 'a2')
R = PolynomialRing(F, 'x', n, order = 'lex') 
x = list(R.gens()) 
a = F.gens()[0]
R1 = PolynomialRing(F, 'x', 3, order = 'lex')
P = (21*a^3 + 73*a^2 + 19*a + 24, 55*a^3 + 48*a^2 + 77*a + 82)
card = E.cardinality()
Q = (30*a^3 + a^2 + 38*a + 99, 56*a^3 + 43*a^2 + 4*a + 81)

(y0, y1) = R2.gens()
f1 = f(y0, y1, xR)
list_1 = weil_descent_1(f1, 2, F)
# print(list_1)
I = R2.ideal(list_1)
print(I)
solution = solve_grobner_basis(I, R2)













#f0 = x^2*y^2 + ( 40)*x^2*y + (15)*x^2 + (40)*x*y^2 + (-36)*x*y + (-23)*x + (15)*y^2 + (-23)*y + (-12)
#f1 = (- 46)*x^2*y + (35)*x^2 + (- 46)*x*y^2 + (31)*x*y + (37)*x + (35)*y^2 + (37)*y + (9)
#f2 = (- 32)*x^2*y + (- 21)*x^2 + (- 32)*x*y^2 + (42)*x*y + (- 5)*x + (- 21)*y^2 + (- 5)*y + (37)
#f3 = (13)*x^2*y + (-21)*x^2 + (13)*x*y^2 + (42)*x*y + (-39)*x + (-21)*y^2 + (-39)*y + (26)
#print([f0, f1, f2, f3])










