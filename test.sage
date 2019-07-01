import time
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
    n = len(R.gens())
    x = R.gens()
    t = list(x)
    solution = []
    print time.time()
    F = R.base_ring()
    I = R.ideal(S)
    basis = I.groebner_basis()
    print time.time()
    elim_ideal = []
    elim_ideal += [elimination_ideal(basis, R, n-1)]
    R1.<X> = F[]
    t[-1] = X
    f = elim_ideal[-1][0](t)
    roots = [[y[0]] for y in f.roots()]
    solution += roots
    print(solution)
    for i in range(n-1)[::-1]:
        J = elimination_ideal(basis, R, i)
        solution1 = []
        for k in J:
            solution2 = solution    
            for j in solution2:
                t1 = list(x)
                t1[i] = X 
                t1[i+1:] = j 
                f = k(t1)
                if f.degree(X) == 0:
                    break
                else:    
                    roots = [[y[0]] for y in f.roots()] 
                    for s in roots:
                        solution1 += [s + j] 
                    solution2.remove(j)
        solution = solution1
        print(solution)            
    print time.time()
    return solution

#F = F_{q^n2}, q = p^n1.    
n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
p = input("insert the characteristic of the base field.")
m = input("number of variables in the polynomial ring.")
Fn = GF(p^n1, 'a')
F, f = Fn.extension(n2, 'a2', map = True)
a = Fn.gens()[0]
a2 = F.gens()[0]
R = PolynomialRing(F, 'x', m, order = 'lex')
(x0, x1, x2, x3) = R.gens()
S = [3*x0**2 + 2*x1*x2 - 2*x0*x3, 2*x0*x2 - 2*x1*x3, 2*x0**2 - 2*x2 - 2*x2*x3, x0**2 + x1**2 + x2**2 - 1]
#input("system of multivariate polynomials to be solved.")
solution = solve_groebner_basis(S, R)
#S = [3*x0**2 + 2*x1*x2 - 2*x0*x3, 2*x0*x2 - 2*x1*x3, 2*x0**2 - 2*x2 - 2*x2*x3, x0**2 + x1**2 + x2**2 - 1]