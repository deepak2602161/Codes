# working in the prime finite field.
p = input("Enter the characteristic of the field.")
F = GF(p)
A, B = input("Enter the parameters of the elliptic curve, separated by comma.")
E = EllipticCurve(F, [A, B])
card = E.cardinality()
r_points = E.rational_points()
P = E.gens()[0]
R = 2*P
l = line([[P[0], P[1]], [R[0], R[1]]], rgbcolor = 'green')
t = text("1", (P[0], P[1]), horizontal_alignment = "left")
for i in range(2, card-1):
    Q1 = i*P
    Q2 = (i+1)*P
    l = l + line([[Q1[0], Q1[1]], [Q2[0], Q2[1]]])
    t = t + text('%s'%i,  (Q1[0], Q1[1]), horizontal_alignment = "left", color = 'black')
M = E.plot(rgbcolor = 'black') + l + t + line([[0, (p)/2], [p-1, p/2]], rgbcolor = 'red')
M.save('tracing_the_generator.png')
print(P.order())