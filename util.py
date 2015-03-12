import math
import random 
from fractions import gcd


def base_coeffs(b,N):
	""" Returns coefficients of N to the base b """
	coeffs = list()

	while N != 0:
		coeffs.append(N % b)
		N /= b

	return coeffs

def is_probable_prime(n):
    """ Miller-Rabin primality test """

    assert n >= 2
    # special case 2
    if n == 2:
        return True
    # ensure n is odd
    if n % 2 == 0:
        return False
    # write n-1 as 2**s * d
    # repeatedly try to divide n-1 by 2
    s = 0
    d = n-1
    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2**s * d == n-1)

    _mrpt_num_trials = 10 
 
    # test the base a to see whether it is a witness for the compositeness of n
    def try_composite(a):
        if pow(a, d, n) == 1:
            return False
        for i in range(s):
            if pow(a, 2**i * d, n) == n-1:
                return False
        return True # n is definitely composite
 
    for i in range(_mrpt_num_trials):
        a = random.randrange(2, n)
        if try_composite(a):
            return False
 
    return True # no base tested showed n as composite

def primes_less_than(N):
    """ Find all primes less than N """ 
        
    sieve = [True] * N

    for i in xrange(3,int(N**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((N-i*i-1)/(2*i)+1)

    return [2] + [i for i in xrange(3,N,2) if sieve[i]]

def is_irreducible(f):
    """ Checks that f is irreducible over the integers """
    pass 


def roots_mod_p(f,p):
    """ Find all roots of a monic irreducible polynomial using Chiens Search  """ 

    pass 

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
    """ Modular inverse of a """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m


	
