# This code is written for zero dimensional ideal.
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

def check_system(S, R):
    n = len(R.gens())
    x = R.gens()
    t = list(x)
    F = R.base_ring()
    I = R.ideal(S)
    G = I.groebner_basis()
    if 1 in G:
        print("corresponding variety is empty.")
        return 0
    for i in t:
        count = 0
        for j in G:
            if (j.lm().is_univariate()) and (j.lm().degree(i) >= 1):
                print(j)
                count += 1
        print(i, count)
        if count == 0:
            print("system is underdetermined and hence the ideal is not zero dimensional(infinite solution in case of infinite field).")
            return 1
    return R.ideal(S)

def solve_groebner_basis(S, R):
    solution = []
    n = len(R.gens())
    x = R.gens()
    t = list(x)
    F = R.base_ring()
    I = check_system(S, R)
    if I == 1:
        print("system is underdetermined and hence the ideal is not zero dimensional(infinite solution in case of infinite field).")
        return 1
    elif I == 0:
        print("corresponding variety is empty.")
        return 0    
    basis = I.groebner_basis()
    # print(basis)
    J = elimination_ideal(basis, R, n-1)
    R1.<X> = F[]
    t[-1] = X
    f = J[-1](t)    
    roots = [[y[0]] for y in f.roots()]
    solution += roots
    # print(solution)
    for i in range(n-1)[::-1]:
        if solution == []:
            return solution
        J = elimination_ideal(basis, R, i)
        # print(J)
        solution1 = []
        for j in solution:
            l = []
            for k in J:
                t1 = list(x)
                t1[i] = X 
                t1[i+1:] = j 
                l += [k(t1)]
            # print(l)
            I1 = R1.ideal(l)
            g = I1.groebner_basis()[0]
            # print(g)        
            roots = [[y[0]] for y in g.roots()]
            for s in roots:
                solution1 += [s + j]
                # print(solution1)
        solution = solution1
        # print(solution)
    return solution  

#F = F_{q^n2}, q = p^n1.    
n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
p = input("insert the characteristic of the base field.")
m = input("number of variables in the polynomial ring.")
if p == 0:
    field = input("rational field or complex field(write either rational or complex).")
    if field == "rational":
        if (n2 > 1) and (n1 > 1):
            modulus = input("insert the modulus of the base number field.")
            modulus1 = input("insert the modulus of the extended number field.")
            Fn = QQ.extension(modulus, 'a')
            F = Fn.extension(modulus1, 'a1')
        elif (n1 > 1) and (n2 == 1):    
            modulus1 = input("insert the modulus of the extended number field.")
            Fn = QQ
            F = Fn.extension(modulus1, 'a1')
        elif (n1 == 1) and (n2 > 1):
            modulus = input("insert the modulus of the base number field.")
            Fn = Numberfield(modulus, 'a')
            F = Fn    
        else:
            Fn = QQ 
            F = QQ   
    else:
        F = CC
else:
    Fn = GF(p^n1, 'a')
    F, f = Fn.extension(n2, 'a2', map = True)
    a = Fn.gens()[0]
    a2 = F.gens()[0]

R = PolynomialRing(F, 'x', m, order = 'lex')
# print(R.gens())
x = list(R.gens())  
S = [x[1] - x[0] - 1, x[1] - x[0] - 1]
# [3*x0**2 + 2*x1*x2 - 2*x0*x3, 2*x0*x2 - 2*x1*x3, 2*x0**2 - 2*x2 - 2*x2*x3, x0**2 + x1**2 + x2**2 - 1]
#input("system of multivariate polynomials to be solved.")
solution = solve_groebner_basis(S, R)
print(len(solution))
print(solution)
# sol = []
# for s in solution:
#     count = 0
#     for f in S:
#         if f(s) == 0:
#             count += 1
#     if count == len(S):
#         sol += [s]
# print(len(sol))           
#S = [3*x0**2 + 2*x1*x2 - 2*x0*x3, 2*x0*x2 - 2*x1*x3, 2*x0**2 - 2*x2 - 2*x2*x3, x0**2 + x1**2 + x2**2 - 1]