import math
import random
from fractions import gcd
import sympy as sp
import numpy as np 
import theano as th 
import pandas as pa 


""" The general number field sieve """
	
def base_coeffs(base,N):
	""" The coefficients the base base expansion of N """
	coefficients = list()
	while N !=0:
		coefficients.append(N%base)
		N/=base 
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
    number_of_trials = 20

    # test the base a for witnesses for the compositeness of N
    def try_composite(a):
        if a ** d % N == 1:
            return False
        for i in range(s):
            if a ** (2 ** i * d) % N == N - 1:
                return False
        return True  # N is definitely composite

    for i in range(number_of_trials):
        a = random.randrange(2, N)
        if try_composite(a):
            return False

    return True  # no base tested showed N as composite

def extended_gcd(aa, bb):
    """ Extended euclidean algorithm """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(
            lastremainder, remainder)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)

def modinv(a, m):
    """ Inverse of a modulo m """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m

def lazy_factors(n):
    """ ALL numbers (not just primes) dividing N """    
    return list(reduce(list.__add__,([i, n//i] \
    	for i in range(1, int(n**0.5) + 1) if n % i == 0)))

def irreducible(f):
    """ Rational roots test of f over Q """
    if f.eval(0) == 0:
        return False 

 	# Leading and constant terms of f
    leading	= [1,-1] + lazy_factors(f.coeffs()[0]) \
        + [-x for x in lazy_factors(f.coeffs()[0])]
    constant= lazy_factors(f.coeffs()[len(f.coeffs()) -1]) \
        + [-x for x in lazy_factors(f.coeffs()[len(f.coeffs()) -1])]

    for q in leading:
        for p in constant:
            if f.eval(p/q) == 0:
                return False 

    return True 

def modInv(a, m):
    """ Modular inverse """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m

def lenstra_elliptic_curve_factor(N):
    """ Lenstra's elliptic curve factoring method """

    if N < 0:
        return [-1] + lenstra_elliptic_curve_factor(-N)
    if N == 0:
        return [0]

    if N == 1:
        return [1] 
    if N == 2 or is_probable_prime(N):
        return [N]

    # Initialize two random integers mod N 
    x0, y0 = random.randrange(1, N), random.randrange(1, N)

    # List of factors to be returned 
    factors = list() 

    # bound for number of trials 
    bound = int(math.sqrt(N))

    for a in xrange(2,N):
        # Build curve out of random points
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
                    return lenstra_elliptic_curve_factor(d) + \
                    lenstra_elliptic_curve_factor(int(N/d))
                else:
                    s = (y - y0) * modInv(x - x0,N)
                    x = s**2 - x - x0  
                    y = - y + s * (s**2 - x - x0 - x)


def legendre(a,p):
    """ Computes the legendre sybol of a mod p """

    # p must be an odd prime 
    if not is_probable_prime(p) or p <=2:
        raise Exception("p must be prime")

    # result 
    s = 1

    def f(q,p):
        """ Cqlculqtes legendre symbol where q is prime """

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
    for x in lenstra_elliptic_curve_factor(a):
        s *= f(x,p)

    return s 


def gnfs(N):
    """ Factors N into list of primes with repitition """

    if is_probable_prime(N):
        return [N]

    # bitlength of N
    k = math.log(N, 2)
    if k >= 110:
        d = 5.0
    elif 80 < k < 110:
        d = 4.0
    else:
        d = 3.0
    
    c = 0
    while True:
        m=int(math.floor(N**(1/d)))+c
        coefficients=base_coeffs(m,N)

        # Create polynomial f
        x = sp.symbols('x')
        f = sp.poly(x, domain=sp.ZZ) - x
        for i in xrange(0, len(coefficients)):
            f += coefficients[i] * x ** i

        # Test if it's irreducible over Q
        if irreducible(f):
            break 
        c +=1 

    # Factor bases for sieving
    rational,algebraic,quadratic = list(),list(),list()

    rational = primes_less_than(m)
    for p in primes_less_than(3 * m):
        for r in sp.ground_roots(f, modulus=p):
            algebraic.append((r, p))
    for p in xrange(3 * m, 4 * m):
        if is_probable_prime(p):
            for r in sp.ground_roots(f, modulus=p):
                quadratic.append((r, p))

    factor_base_length = len(rational)+len(algebraic)+len(quadratic)
    U = np.zeros(shape=(factor_base_length+2,factor_base_length+1))

    number_of_smooths = 0
    sample_size = 10000 
    while number_of_smooths < factor_base_length:
        b = random.randrange(1, m)

        S = xrange(0,sample_size)
        R = [[] for _ in S]
        A = [[] for _ in S]

        # find smooth integers in Z and Z[theta]
        for i in xrange(0,len(S)):
            is_rational_smooth  = False  
            is_algebraic_smooth = False 

            for p in rational:
                if (S[i] + b*m) %p == 0:
                    l=1
                    while (S[i] + b*m) %(p**(l+1)) == 0: 
                        l+=1
                    R[i].append((p,l))

            prod = 1 
            for (p,l) in R[i]:
                prod*=p**l 
            if S[i]+b*m == prod:
                is_rational_smooth = True 

            for (r,p) in algebraic:
                if (S[i] + b*r) %p == 0:
                    l=1 
                    # while (S[i] + b*r) %(p**(l+1)) == 0:
                    #     l+=1
                    A[i].append(((r,p),l))

            prod = 1 
            for ((r,p),l) in A[i]:
                while l > 0:
                    prod*=p
                    l-=1
            if prod == (-1)**(sp.degree(f))*f.eval(-S[i]/b):
                is_algebraic_smooth = True 

            if is_rational_smooth and is_algebraic_smooth:
                for j in xrange(0,len(rational)):
                    p = rational[j]
                    l=0
                    while (S[i] + b*m) %(p**(l+1)) == 0: 
                        l+=1
                    U[number_of_smooths,j] = l % 2

                for j in xrange(0,len(algebraic)):
                    (r,p) = algebraic[j]
                    l=0
                    if (S[i]+b*r) % p == 0:
                        l=1
                    # while (S[i]+b*r) %(p**(l+1)) == 0:
                        # l+=1  
                    U[number_of_smooths,len(rational)+j] =  l % 2

                for j in xrange(0,len(quadratic)):
                    (s,q) = quadratic[j]
                    if legendre(S[i]+b*s,q) == 1:
                        U[number_of_smooths,len(rational)+len(algebraic)] = 0
                    else:
                        U[number_of_smooths,len(rational)+len(algebraic)] = 1 
                    
                number_of_smooths +=1
                
    np.Solve(U)

