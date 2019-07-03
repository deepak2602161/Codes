# Codes
Codes for ECDLP
I have written the above codes to solve the ECDLP problem using Index Calculus Algorithm given by Gaudry(Index calculus for abelian varieties of small dimension and the elliptic curve discrete logarithm problem).
1. addition_elliptic_curve_geometric.src
This code was written to see how addition operation in the elliptic curve group works geometrically.
Input Format: characteristic of the field // parameters of the elliptic curve separated by comma.
Output : a .png file shows the line starting from the point P and ends at the reflection of the output of the points P and Q.

2. tracing_the_generator.src
This code was motivated by the above code and to understand the complexity of the group operation using geometry.
Input Format: characteristic of the field // parameters of the elliptic curve separated by comma.
Output: a .png file showing how a generator traverses through the group geometrically.

3. better_brute_force.src:
This code was motivated by the simple formula: o(q) = n/gcd(n, j), where G is a cyclic group with a generator g such that 
o(g) = n, and q is any arbitrary element and j is the discrete log base g. This idea behind the code was to reduce the search 
space for the brute force in case q is not a generator, which can be a case if the n is not prime. This code will work only if 
the element q is not the generator, otherwise there is not point of this code.

4. line_finite_field.src:
This code contains a function which has arguments (N, m, c, P, R).
N: size of the characteristic of the finite field,
m: slope of the line,
c: y-intercept of the line,
P: starting point of the line,
R. end point of the line.

5. summation_polynomial.src:
Input: field, parameters of the elliptic curve, and the index(k) for the summation polynomial, k>2.
Ouput: Sk

6. zero_multi_system.sage:
Input format: this code will take the field, and the system of 'n' multivariate polynomial equations in n variables.
Output: corresponding affine variety of the above given ideal.
This code is written to solve system of multivariate polynomial equations for the zero dimensional ideal only, but works for 
any field.

7. solve_system_multi_equation.sage:
This code will work for finite field only but for any ideal. Word is under progress.

8. dlog_gaudry1.src:
This code is written to solve the discrete log problem as promised above. Work is under progress.
