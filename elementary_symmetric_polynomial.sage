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
    
n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
p = input("insert the characteristic of the base field.")
n = input("number of variables in the polynomial ring.")
Fn = GF(p^n2, 'a')
F = Fn.extension(n2, 'a2')
R = PolynomialRing(F, 'x', n, order = 'lex') 
e = elementary_symmetric_polynomials(n, R)
print(e)    