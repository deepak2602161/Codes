import time
#k should be greater than 2.
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

n = input("insert the degree of the extension field")
p = input("insert the prime.")
k = input("insert the index for the required summation polynomial.")
A, B = input("insert the parameters of required elliptic curve.")
F = GF(p^n, 'a')
R = PolynomialRing(F, 'x', k, order = 'lex')
x = list(R.gens())
print time.time()
print(summation_poly(k, R, A, B).degree(x[0]))
print(summation_poly(k, R, A, B))
print time.time()