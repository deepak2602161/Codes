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

n1 = input("insert the degree of the extension field.")
n2 = input("insert the degree of the extension of the base field.")
p = input("insert the characteristic of the base field.")
n = input("number of variables in the polynomial ring.")
Fn = GF(p^n2, 'a')
F = Fn.extension(n2, 'a2')
R = PolynomialRing(F, 'x', n, order = 'lex')
rep = input("insert the symmetric multivariate polynomial you want to convert. input format = [[a_11, a_12, ..., a_1n, c_1], ..., [a_n1, a_n2, ..., a_nn, c_n]], so f(X) = c_1*x^alpha_1, ..., c_n*x^alpha_n, where alpha_i = (a_i1, ..., a_in)")     
f = dict_to_poly(rep, R)
g = conversion_to_esp(f, R)
print(g)