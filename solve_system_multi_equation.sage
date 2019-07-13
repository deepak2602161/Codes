# This code is written for SMVP over finite field.
import time
from sage.rings.polynomial.polydict import PolyDict
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

#F = F_{q^n2}, q = p^n1(Finite Field) or number field or complex field.    
n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
p = input("insert the characteristic of the base field.")
n = input("number of variables in the polynomial ring.")
Fn = GF(p^n2, 'a')
F = Fn.extension(n2, 'a2')
R = PolynomialRing(F, 'x', n, order = 'lex') 
x = list(R.gens()) 
S = input("insert system of multivariate polynomials to be solved. input format = [[[a_11, a_12, ..., a_1n, c_1], ..., [a_n1, a_n2, ..., a_nn, c_n]], ...,], so f_1(X) = c_1*x^alpha_1, ..., c_n*x^alpha_n, where alpha_i = (a_i1, ..., a_in)")   
S = [f, g]
print time.time()
solution = solve_groebner_basis(S, R)
print time.time()
print(len(solution))
print(solution)