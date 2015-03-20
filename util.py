#!/usr/bin/env python

# TODO: 
# - Fix bugs in factor
# - Test Chiens vs Factoring method 
# - Miller rabin irreduciblity 
# - Fix lenstra, 1 million or less... 

import math
import random 
from fractions import gcd
from sympy import * 

""" 
Many methods for the general number field sieve 

"""

def extended_gcd(aa, bb):
    """ Extended Euclidean algorithm """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0

    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
        y, lasty = lasty - quotient*y, y

    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)
 

def modinv(a, m):
    """ Modular inverse of a modulo m """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError

    return x % m

def is_probable_prime(N):
    """ Miller-Rabin primality test on an integer N """
    assert N >= 2
    # special case 2
    if N == 2:
        return True
    # ensure N is odd
    if N % 2 == 0:
        return False
    # write N-1 as 2**s * d
    # repeatedly try to divide N-1 by 2
    s = 0
    d = N-1
    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2**s * d == N-1)

    # This needs to be set, but 10 corresponds to false possitive probablility = 1/2**10
    _mrpt_num_trials = 10 
 
    # test the base a to see whether it is a witness for the compositeness of N
    def try_composite(a):
        if a**d % N == 1:
            return False
        for i in range(s):
            if a**(2**i*d) % N == N-1:
                return False
        return True # N is definitely composite
 
    for i in range(_mrpt_num_trials):
        a = random.randrange(2, N)
        if try_composite(a):
            return False
 
    return True # no base tested showed N as composite


def lenstra_elliptic_curve_factor(N):
    """ Lenstra's elliptic curve factoring method """

    # Base case 
    if N == 1:
        return []

    # Cant factor a prime! 
    if N == 2 or is_probable_prime(N):
        return [N]

    # Initialize two random integers mod N 
    x0, y0 = random.randrange(1, N), random.randrange(1, N)
    factors = list() 
    bound = int(math.sqrt(N))

    for a in xrange(2,N):
        # Build curve our of random points
        b = y0**2 - x0**3 - a*x0

        # Check curve is not singular 
        if 4*a**3 - 27*b**2 ==0:
            next 

        # Initially double point 
        s = (3*x0**2 + a) 
        (x,y) = (s**2 - 2*x0, s*((s**2 - 2*x0) - x0) - y0)

        # Keep adding points until gcd(x-x0,N) != 1
        for k in xrange(2,bound):
            for i in xrange(0,math.factorial(k)):
                d = gcd(x- x0,N)
                if d != 1:
                    return [d] + lenstra_elliptic_curve_factor(N/d)
                else:
                    s = (y - y0) * modinv(x - x0,N)
                    x = s**2 - x - x0  
                    y = - y + s * (s**2 - x - x0 - x)


def base_coeffs(b,N):
	""" Returns coefficients of N to the base b """
	coeffs = list()

	while N != 0:
		coeffs.append(N % b)
		N /= b

	return coeffs


def is_probable_prime(N):
    """ Miller-Rabin primality test on an integer N """
    assert N >= 2
    # special case 2
    if N == 2:
        return True
    # ensure N is odd
    if N % 2 == 0:
        return False
    # write N-1 as 2**s * d
    # repeatedly try to divide N-1 by 2
    s = 0
    d = N-1
    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2**s * d == N-1)

    # This needs to be set, but 10 corresponds to false possitive probablility = 1/2**10
    _mrpt_num_trials = 10 
 
    # test the base a to see whether it is a witness for the compositeness of N
    def try_composite(a):
        if pow(a, d, N) == 1:
            return False
        for i in range(s):
            if pow(a, 2**i * d, N) == n-1:
                return False
        return True # N is definitely composite
 
    for i in range(_mrpt_num_trials):
        a = random.randrange(2, N)
        if try_composite(a):
            return False
 
    return True # no base tested showed N as composite


def primes_less_than(N):
    """ Find all primes less than N """ 
    sieve = [True] * N

    for i in xrange(3,int(N**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((N-i*i-1)/(2*i)+1)

    return [2] + [i for i in xrange(3,N,2) if sieve[i]]


def Chien_Search(f,p):
    """ 
    Chien's search for roots of f mod p:
    http://en.wikipedia.org/wiki/Chien_search
    """ 

    # to represent roots in terms of primitive element 
    a = primitive_root(p)
    coefficients = f.coeffs()
    A = coefficients[::-1]
    roots = list()

    if f.eval(0) % p == 0:
        roots.append(0)

    if sum(A) % p ==0: 
        roots.append(1)

    for i in xrange(0,p-1): 
        for j in xrange(1,len(coefficients)):
            A[j] = A[j]*a**i 
        if sum(A) % p == 0:
            roots.append(a**i % p)

    return roots       


def rand_poly(p,n):
    """ Random polynomial with degree n and coefficients in finite field Z_q """

    # Set degree of random polynomial 
    d = random.randrange(1, n)

    # Indeterminant variable x 
    x = symbols('x')

    # zero polynomial in F_p[x]
    f = poly(0,x,modulus=p)

    for i in xrange(0,d):
        a = random.randrange(0, p)
        f = poly(f+ a*x**i,x,modulus=p) 

    return f 


def square_free_factorization(f,p):
    """ The square free factorizaiton of f mod p """

    # Initial variables
    i = 1
    r = 1 
    g = diff(f) # Derivative of f wrt x 

    # Alias square_free_factorization 
    F = square_free_factorization

    if g != 0:
        c = gcd(f,g)
        w = f / c 
        while w != 1: 
            y = gcd(w,c)
            z = w/y 
            r = r*z**i
            i += 1 
            w = y 
            c = c/y 
        if c!=1: 
            c = c**(1/p)
            return r*F(c,p)**p 
        else:
            return r 
    else:
        f = f**(1/p)
        return F(f,p)**p 


def distinct_degree_factorization(f,p):
    """ 
    Finds all pairs (g, d), such that 
        -f has an irreducible factor of degree d and
        -g is the product of all monic irreducible factors of f of degree d. 

    """

    # Initial variables 
    i = 1 
    S = []
    ff = f

    # Indeterminant variable x 
    x = symbols('x')

    while degree(ff) >= 2*i:
        g = gcd(ff,x**(p*i) - x)
        if g != 1:
            S.append((g,i))
            ff /= g 
        i += 1 

    if ff !=1:
        S.append((ff,degree(ff)))
    if S == []:
        return [(f,1)]
    else:
        return S 


def cantor_zassenhaus_algorithm(f,p,d):
    """ 
    Cantor-Zassenhaus's algorithm for finding 
    the set of monic irreducible factors of f
    """

    Factors = [f]
    n = degree(f)
    r = n / d 

    while len(Factors) < r:
        h = rand_poly(p,n - 1)
        g = reduce(h ** ((p**d - 1)/2) - 1, [f])
        for u in [x for x in Factors if degree(x) > d]:
            if gcd(g,u) != 1 and gcd(g,u) != u:
                Factors = Factors.remove(u) + [gcd(gcd(g,u),u/gcd(g,u))]

    return Factors

def factor(f,p):
    """ Factor univariate polynomial in a finite field of order p """

    # list of factors to return
    Factors = list()

    # f must be monic 
    assert f.coeffs()[degree(f)-1] == 1 

    # Square free factorization of f
    f1 = square_free_factorization(f,p)  

    # all pairs (g,d) satisfying 
    # f has an irreducible factor of degree d and
    # g is the product of all monic irreducible factors of f of degree d.
    G = distinct_degree_factorization(f1,p)

    # The set of monic irreducible factors of f
    for (g,d) in G:
        Factors += cantor_zassenhaus_algorithm(f,p,d) 

    return Factors 


def primitive_root(p):
    """ The first primitive root of p """

    # prime factors of p-1 
    phi_factors = lazy_factors(p-1)

    # test if a has order less than p-1 
    for a in xrange(2,p):
        primitive = True  
        for l in phi_factors:
            if a**((p-1)/l) % p == 1:
                primitive = False 
        if primitive:
            return a


def lazy_factors(N):    
    """ Simple factoring meth for small integers """
    return set(reduce(list.__add__, 
                ([i, N//i] for i in range(1, int(N**0.5) + 1) if N % i == 0)))


	
