import time
def non_prime_cyclic_curve(i, A,  B):
    p = next_prime(i)
    F = GF(p)
    E = EllipticCurve(F, [A, B])
    card = E.cardinality()
    while is_prime(card) == 1:
        p = next_prime(p)
        F = GF(p)
        E = EllipticCurve(F, [A, B])
        card = E.cardinality()
    return [p, E]   

def generator(p, E, A, B):
    P = E.gens()[0]
    card = E.cardinality()
    while P.order() != card:
        E = non_prime_cyclic_curve(p, A, B)[1]
        p = non_prime_cyclic_curve(p, A, B)[0]
        card = E.cardinality()
        P = E.gens()[0]
    return [p, E, P]    

Z = IntegerRing()
i = 100000
E = non_prime_cyclic_curve(i, -3, 1)[1]
p = non_prime_cyclic_curve(i, -3, 1)[0]
r = E.cardinality()
print time.time()
[p, E, P] = generator(p, E, -3, 1)
r = E.cardinality()
print time.time()
print(E, P, r, P.order())


a = randrange(r)
Q = a*P
t = Q.order() 
print time.time()
print('a, r, t', a, r, t, Q)

def possible_dlog(t, r):
    d = r/t
    Zn = Integers(r/d)
    l1 = [a for a in Zn if gcd(a, r/d) == 1]
    l = []
    for i in range(len(l1)):
        l = l + [d*Z(l1[i])]
    return l

l = possible_dlog(t, r)
if t != r:    
    for i in range(len(l)):
        if l[i]*P == Q:
            print('dlog', l[i])
            break
print time.time()