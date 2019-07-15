import time
from sage.rings.polynomial.polydict import PolyDict

def elementary_symmetric_polynomials(n, R):
    F = R.base_ring()
    x = list(R.gens())
    e = (n+1)*[0]
    l = monomials(x, n*[2])
    l_degree = [y.degree() for y in l]
    for t in l:
        t_degree = t.degree()
        e[t_degree] += t     
    return e

def conversion_to_esp(f, R):
    n = len(R.gens())
    x = list(R.gens())
    e = elementary_symmetric_polynomials(n, R)
    g = 0
    h = f
    while h != 0:
        m = h.lm()
        m_exp = m.exponents()[0]
        t1 = h.lc()
        t2 = h.lc()   
        for j in range(n-1):
            t1 = t1*(x[j]**(m_exp[j] - m_exp[j+1]))
            t2 = t2*(e[j+1]**(m_exp[j] - m_exp[j+1]))
        t1 = t1*(x[n-1]**(m_exp[n-1]))
        t2 = t2*(e[n]**(m_exp[n-1]))
        g += t1
        h = h - t2  
    return g 

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
def curve(i, n1, n2, A,  B):
    p = next_prime(i)
    Fn = GF(p^n2, 'a')
    F = Fn.extension(n1, 'a1')
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
    Fn = F.base_ring()
    n2 = Fn.degree()
    p = F.characteristic()
    a = Fn.gens()[0]
    factorbase = []
    z = 0
    if F(B).is_square() == 1:
        t = F(B).nth_root(2)
        P1 = (z, t)
        factorbase += [P1]
    
    if n2 == 1:
        for i in range(p^n2):
            z = i*a
            if F(z**3 + A*z + B).is_square() == 1:
                t = F(z**3 + A*z + B).nth_root(2)
                P1 = (z, t)
                factorbase += [P1]
    else:
        for i in range(p^n2):
            z = a**i
            if F(z**3 + A*z + B).is_square() == 1:
                t = F(z**3 + A*z + B).nth_root(2)
                P1 = (z, t)
                factorbase += [P1]           
    return factorbase

# eliptic curve discrete log computation using index calculus algorithm given by Gaudry.
def relations(P, Q, F, A, B):
    count = 0 #number of relations
    factorbase = factor_base(F, A, B)
    s = len(factorbase)
    print(s)
    M = matrix(F, s+1, s)
    a = matrix(F, s+1, 1)
    b = matrix(F, s+1, 1)
    E = EllipticCurve(F, [A, B])
    r = E.cardinality()
    print(r)
    Fn = F.base_ring()
    n2 = Fn.degree()
    n1 = F.degree()
    p = F.characteristic()
    a = Fn.gens()[0]
    R = PolynomialRing(F, 'x', n1+1, order = 'lex')
    x = list(R.gens())
    f = summation_poly(n1 + 1, R, A, B)
    while (count < s+1):
        u = randrange(r)
        v = randrange(r)
        R = u*E(P) + v*E(Q) 
        t = x[:n1] + [R[0]]
        f1 = f(t)
        print(f1)
        f1_esp = conversion_to_esp(f1, R)
        print(f1_esp)
        count = s+1




#F = F_{q^n2}, q = p^n1(Finite Field) or number field or complex field.    
n1 = 4
#input("insert the degree of the extension field.")
n2 = 1
#input("insert the degree of the extension of the base field.")
p = 1009
F = GF(p^n1, 'a')
a = F.gens()[0]
A, B = [529*a^3 + 210*a^2 + 379*a + 351, 636*a^3 + 595*a^2 + 7*a + 216]
#input("insert the size of the prime.")
#input("insert the parameters of required elliptic curve.")
P = (748*a^3 + 600*a^2 + 187*a + 357, 322*a^3 + 347*a^2 + 734*a + 656)
Q = (60*a^3 + 373*a^2 + 429*a + 954, 634*a^3 + 528*a^2 + 537*a + 145)
relations(P, Q, F, A, B)

#p = 1009, n1 = 4, n2 = 1, A = 529*a^3 + 210*a^2 + 379*a + 351, B = 636*a^3 + 595*a^2 + 7*a + 216, 
#P = (748*a^3 + 600*a^2 + 187*a + 357 : 322*a^3 + 347*a^2 + 734*a + 656 : 1)
#Q = (60*a^3 + 373*a^2 + 429*a + 954 : 634*a^3 + 528*a^2 + 537*a + 145 : 1), Q = a*P, a = 504988480210















