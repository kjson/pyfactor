import math
import random 


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
    """ Find all primes less than a given bound """ 
        
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

    # find primiteve root of p 

    phi_factors = easy_factor(p-1)

    for a in xrange(2,p):

        primitive = True  

        for l in phi_factors:
            if a**((p-1)/l) % p == 1:
                primitive = False 

        if primitive:
            return a


def easy_factor(N):

    """ Dirty Lenstras EC factoring algorithm """

    x_0 = random.randrange(2, N)
    y_0 = random.randrange(2, N)

    bound = int(math.sqrt(N))

    def extended_gcd(aa, bb):
        lastremainder, remainder = abs(aa), abs(bb)
        x, lastx, y, lasty = 0, 1, 1, 0
        while remainder:
            lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
            x, lastx = lastx - quotient*x, x
            y, lasty = lasty - quotient*y, y
        return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)

    def modInv(a, m):
        g, x, y = extended_gcd(a, m)
        return g,x % m 

    def add(a,N,(x1,y1),(x2,y2)):

        if x1 != x2:  
            d, x = modInv(x1 - x2,N) 
            if d != 1:
                return d, (1,1)

            s = (y1 - y2) * x

        else:
            s = (3*x1**2 + a) * modInv(2*y1,N) 

        x3 = s**2 - x1 - x2 
        y3 = - y1 + s * (x3 - x1)

        return 1,(x3 % p,y3 % N)

    def scale(a,N,num,(x1,x2)):
        if num % 2 == 0:
            return scale(a,N,num / 2, add(a,N,(x1,x2),(x1,x2)))
        else:
            return add(a,N,(x1,x2), scale(a,N,num - 1,(x1,x2)))
        
    for a in xrange(2,N):
        b = y_0**2 - x_0**3 - a*x_0

        if 4*a**3 - 27*b**2 ==0:
            continue

        for k in xrange(2,bound):
            if scale(a,N,k,(x_0,y_0))[0] != 1:
                if gcd(N,scale(a,N,k,(x_0,y_0))[1]) != 1:
                    return gcd(N,scal(a,N,k,(x_0,y_0))[1])



roots_mod_p(3,32)





	
