import math

from fractions import gcd

import random as rn
import sympy as sp
import numpy as np


""" The general number field sieve. """


def base_coeffs(b, N):
    """ The coefficients in the base b expansion of N """
    coefficients = list()
    while N != 0:
        coefficients.append(N % b)
        N /= b
    return coefficients


def primes_less_than(N):
    """ Primes less than N """
    sieve = [True] * N
    for i in xrange(3, int(N ** 0.5) + 1, 2):
        if sieve[i]:
            sieve[i * i::2 * i] = [False] * ((N - i * i - 1) / (2 * i) + 1)
    return [2] + [i for i in xrange(3, N, 2) if sieve[i]]


def is_probable_prime(N):
    """ Miller-Rabin primality test on an integer N """
    if N == 2:
        return True
    # ensure N is odd
    if N % 2 == 0:
        return False
    # write N-1 as 2**s * d
    # repeatedly try to divide N-1 by 2
    s = 0
    d = N - 1

    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2 ** s * d == N - 1)

    # hard variable, 10 means false possitive prob = 1/2^{10}
    number_of_trials = 10

    # test the base a for witnesses for the compositeness of N
    def try_composite(a):
        if a ** d % N == 1:
            return False
        for i in range(s):
            if a ** (2 ** i * d) % N == N - 1:
                return False
        return True  # N is definitely composite

    for i in range(number_of_trials):
        a = rn.randrange(2, N)
        if try_composite(a):
            return False

    # no base tested showed N as composite
    return True


def extended_gcd(a, b):
    """ Extended euclidean algorithm """
    lastremainder, remainder = abs(a), abs(b)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(
            lastremainder, remainder)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return lastremainder, lastx * (-1 if a < 0 else 1), \
        lasty * (-1 if b < 0 else 1)


def mod_inv(a, m):
    """ Inverse of a modulo m """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m


def lazy_factors(n):
    """ ALL numbers (not just primes) dividing N """
    return list(reduce(list.__add__, ([i, n//i]
                for i in range(1, int(n**0.5) + 1) if n % i == 0)))


def irreducible(f):
    """ Rational roots test of f over Q """
    if f.eval(0) == 0:
        return False

    # Leading and constant terms of f
    leading = [1, -1]+lazy_factors(f.coeffs()[0]) \
        + [-x for x in lazy_factors(f.coeffs()[0])]
    constant = lazy_factors(f.coeffs()[len(f.coeffs()) - 1]) \
        + [-x for x in lazy_factors(f.coeffs()[len(f.coeffs()) - 1])]

    # Check of roots 
    for q in leading:
        for p in constant:
            if f.eval(p/q) == 0:
                return False

    return True

def legendre(a, p):
    """ Computes the legendre symbol of a mod p """
    # p must be an odd prime 
    if not is_probable_prime(p) or p <=2:
        raise Exception("p must be prime")

    # result 
    s = 1

    def f(q,p):
        """ Calculates legendre symbol where q is prime """

        # qll the bqse cqses 
        if q == 1:
            return 1 
        if q % p ==0:
            return 0 
        if q == -1:
            if p % 4 == 1:
                return 1
            else:
                return -1
        if q == 2:
            if p % 8 == 1 or p % 8 == 7:
                return 1
            if p % 8 == 3 or p % 8 == 5:
                return -1

        # Recursive step 
        if q > p:
            return f(q % p,p)

        if q % 4 == 1 or p % 4 == 1:
            return f(p,q)
        else:
            return -f(p,q) 

    # Because the legendre symbol is multiplicative 
    for x in lenstra(a):
        s *= f(x, p)

    return s 

def lenstra(N):
    """ Lenstra's elliptic curve factorization method """

    # If N is negative 
    if N < 0:
        return [-1] + lenstra(-N)

    # If N is prime 
    if 0 <= N <= 3 or is_probable_prime(N):
        return [N]

    # Initialize two random points mod N 
    x0,y0 = rn.randrange(1,N),rn.randrange(1,N)

    # Bound for the number of trials 
    bound = int(math.sqrt(N))

    for a in xrange(1,N):
        # build elliptic curve using random points 
        b = y0**2 - x0**3 - a*x0 

        # Check if the curve is singular 
        if 4*a**3 - 27*b**2 ==0:
            next 

        d = gcd(2*y0,N)
        if d!= 1:
            return lenstra(d) + lenstra(N/d)

        # Initial double of the point (x0,y0)
        s = (3*x0**2 + a) * mod_inv(2*y0,N)
        x,y = (s**2 - 2*x0,s*(3*x0 - s**2)-y0)

        # Search for non-trivial gcd's 
        for k in xrange(1,bound):
            for i in xrange(1,k):
                d = gcd(x-x0,N)
                if d != 1:
                    return lenstra(d) + lenstra(N/d)
                else:
                    s = (y-y0) * mod_inv(x-x0,N)
                    y = (s*(2*x + x0 - s**2)-y) % N
                    x = (s**2 - x -x0) % N


def norm(f,theta,alpha):
    """ Defines a norm on Z[theta]"""
    roots = sp.solve(f)
    norm = 1

    for theta_i in roots:
        norm*= alpha.subs(theta,theta_i)

    return norm.simplify()

def order_p(N,p):
    """ The exponent of p in the prime factorization of N """
    if N % p != 0:
        return 0 
    else:
        e = 1 
        while N % p**(e+1) == 0:
            e+=1
    return e  
        

def gnfs(N):
    """ Factors N into list of primes with repitition """

    # Can't factor a prime number! 
    if is_probable_prime(N):
        return [N]

    # Determine the degree of the polynomial to be used 
    # based on the bitlength of N
    k = math.log(N, 2)
    if k >= 110:
        d = 5.0
    elif 80 < k < 110:
        d = 4.0
    else:
        d = 3.0

    # Select a monic irreducible polynomial of degree d
    # with integer coefficients.  
    c = 0
    while True:
        # Find the coefficients in the m base expansion of N
        # where m is the dth root of N, plus some constant c  
        m=int(math.floor(N**(1/d)))+c
        coefficients=base_coeffs(m,N)

        # Create polynomial f of degree d with coefficients 
        # from 'coefficients'
        x = sp.symbols('x')
        f = sp.poly(x, domain=sp.ZZ) - x
        for i in xrange(0, len(coefficients)):
            f += coefficients[i] * x ** i

        # Test if f is irreducible over Q and that m is 
        # still a root mod N 
        if irreducible(f) and f.eval(m) % N == 0:
            break 
        else:
            c +=1

    # Create rational, algebraic and quadratic factor  
    # bases for sieving
    rational,algebraic,quadratic = list(),list(),list()
    rational = primes_less_than(m)
    for p in primes_less_than(3 * m):
        for r in sp.ground_roots(f, modulus=p):
            algebraic.append((r, p))
    for p in xrange(3 * m, 4 * m):
        if is_probable_prime(p):
            for r in sp.ground_roots(f, modulus=p):
                quadratic.append((r, p))    

    # Create a numpy matrix to store "smooth" elements
    factor_base_length = len(rational)+len(algebraic)+len(quadratic)
    U = np.zeros(shape=(factor_base_length+2,factor_base_length+1))


    # When number_of_smooths > factor_base_length we
    # stop sieving for "smooth" elements
    number_of_smooths = 0
    sample_size = 10000

    # Sieve for smooth pairs (a,b)
    b = 0
    while number_of_smooths < factor_base_length +1:
        # b = rn.randrange(1, m)
        b+=1

        S = xrange(0,sample_size)
        R = [[] for _ in S]
        A = [[] for _ in S]

        # find smooth integers in Z and Z[theta]
        for i in xrange(1,len(S)):
            a = S[i]

            # a and b must be coprime 
            if gcd(a,b) != 1:
                next 

            # Two flags to keep track of whether (a,b) is smoooth in 
            # Z and in Z[theta]
            is_rational_smooth  = False  
            is_algebraic_smooth = False 

            # Build list of primes dividing a + b*m 
            for p in rational:
                if (a + b * m) % p == 0:
                    l=1
                    while (a + b*m) % (p**(l+1)) == 0: 
                        l+=1
                    R[i].append((p,l))

            # Build list of prime ideals "dividing" a + b*theta 
            for (r,p) in algebraic:
                if (a + b*r) % p == 0:
                    l=1 
                    # while (a + b*r) % (p**(l+1)) == 0:
                    #     l+=1 
                    A[i].append(((r,p),l))

            # Check if a + b*m is smooth over the rational factor base 
            # and if a + b*theta is "smooth" over the algebraic factor base 
            prod = 1 
            for (p,l) in R[i]:
                prod*=p**l 
            if a+b*m == prod:
                is_rational_smooth = True 

            prod = 1 
            for ((r,p),l) in A[i]:
                while l > 0:
                    prod*=p
                    l-=1
            if prod == (-1)**(sp.degree(f))*f.eval(-a/b):
                is_algebraic_smooth = True 

            
            # If a+b*m and a+b*theta are smooth over the 
            # rational and algebraic factor bases respectively 
            # then add then add the exponents of their factorizations 
            # to as a row in U.  
            if is_rational_smooth and is_algebraic_smooth:

                # Exponents for rational factor base mod 2 
                for j in xrange(0,len(rational)):
                    p = rational[j]
                    l=0
                    while (a + b*m) %(p**(l+1)) == 0: 
                        l+=1
                    U[number_of_smooths,j] = l % 2

                # Expoenents of algebraic facctor base mod 2 
                for j in xrange(0,len(algebraic)):
                    (r,p) = algebraic[j]
                    l=0
                    if (a+b*r) % p == 0:
                        l=1
                        # while (a+b*r) %(p**(l+1)) == 0:
                        #     l+=1  
                    U[number_of_smooths,len(rational)+j] =  l % 2

                # Quadratic residues
                for j in xrange(0,len(quadratic)):
                    (s,q) = quadratic[j]
                    if legendre(a+b*s,q) == 1:
                        U[number_of_smooths,len(rational)+len(algebraic)] = 0
                    else:
                        U[number_of_smooths,len(rational)+len(algebraic)] = 1 
                    
                number_of_smooths +=1


