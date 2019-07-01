def non_prime_curve(i, A,  B):
    p = next_prime(i)
    F = GF(p)
    E = EllipticCurve(F, [A, B])
    card = E.cardinality()
    while is_prime(card) != 1:
        p = next_prime(p)
        F = GF(p)
        E = EllipticCurve(F, [A, B])
        card = E.cardinality()
    return [p, E]    