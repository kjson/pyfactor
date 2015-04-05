import math
from sympy import *
from fractions import gcd
import random
import numpy as np

""" The number field sieve """

def base_coeffs(b, N):
    """ Returns coefficients of N to the base b """
    coeffs = list()
    while N != 0:
        coeffs.append(N % b)
        N /= b
    return coeffs


def primes_less_than(N):
    """ Find all primes less than N """
    sieve = [True] * N
    for i in xrange(3, int(N ** 0.5) + 1, 2):
        if sieve[i]:
            sieve[i * i::2 * i] = [False] * ((N - i * i - 1) / (2 * i) + 1)
    return [2] + [i for i in xrange(3, N, 2) if sieve[i]]


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
    d = N - 1
    while True:
        quotient, remainder = divmod(d, 2)
        if remainder == 1:
            break
        s += 1
        d = quotient
    assert(2 ** s * d == N - 1)

    # This needs to be set, but 10 corresponds to false possitive probablility
    # = 1/2**10
    _mrpt_num_trials = 10

    # test the base a to see whether it is a witness for the compositeness of N
    def try_composite(a):
        if a ** d % N == 1:
            return False
        for i in range(s):
            if a ** (2 ** i * d) % N == N - 1:
                return False
        return True  # N is definitely composite

    for i in range(_mrpt_num_trials):
        a = random.randrange(2, N)
        if try_composite(a):
            return False

    return True  # no base tested showed N as composite


def extended_gcd(aa, bb):
    """ Extended Euclidean algorithm """
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(
            lastremainder, remainder)
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y
    return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)


def modinv(a, m):
    """ Modular inverse of a modulo m """
    g, x, y = extended_gcd(a, m)
    if g != 1:
        raise ValueError
    return x % m


def lazy_factors(n):
    """ Returns list the is ALL numbers (not just primes) dividing N """    
    return list(reduce(list.__add__,([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))


def irreducible(f):
    """ Test irreducibility of f over Q """

    if f.eval(0) == 0:
        return False 

    leading     = [1,-1] + lazy_factors(f.coeffs()[0]) + [-x for x in lazy_factors(f.coeffs()[0])]
    constant    = lazy_factors(f.coeffs()[len(f.coeffs()) -1]) + [-x for x in lazy_factors(f.coeffs()[len(f.coeffs()) -1])]

    for q in leading:
        for p in constant:
            if f.eval(p/q) == 0:
                return False 

    return True 


# def nfs(N):
#     """ Factors N into list of primes with repitition 
#     TODO:
#     -- Sieve and fill in the array.
#     """

#     # Don't waste time
#     if is_probable_prime(N):
#         return [N]

#     # bitlength of N
#     k = math.log(N, 2)
#     if k >= 110:
#         d = 5.0
#     elif 80 < k < 110:
#         d = 4.0
#     else:
#         d = 3.0

#     while True:
#         i = 0
#         # Calc m such that m^d = N
#         m = int(math.floor(N ** (1 / d))) + i

#         # Coefficients in base-m rep of N
#         coeffs = base_coeffs(m, N)

#         # Create polynomial
#         x = symbols('x')
#         f = poly(x, domain=ZZ) - x
#         for i in xrange(0, len(coeffs)):
#             f += coeffs[i] * x ** i

#         if irreducible(f):
#             break

#     # Factor bases for sieving
#     rational = primes_less_than(m)
#     algebraic = list()
#     quadratic = list()

#     for p in primes_less_than(3 * m):
#         for r in ground_roots(f, modulus=p):
#             algebraic.append((r, p))
#     for p in xrange(3 * m, 4 * m):
#         if is_probable_prime(p):
#             for r in ground_roots(f, modulus=p):
#                 quadratic.append((r, p))

#     # Number of entries to check
#     M = 1000

#     # Random starting value
#     b = random.randrange(1, m)

#     # List of M "soil" integers
#     S = [(a, b) for a in xrange(0, M)]

#     # Lists to keep track of smooth pairs
#     R = [[] for _ in S]
#     A = [[] for _ in S]

#     # Matrix of all smooth pairs
#     U = np.zeros(shape=(len(rational)+len(algebraic)+len(quadratic)+2,len(rational)+len(algebraic)+len(quadratic)+1))

#     # Sieve with rational factor base
#     for p in rational:
#         i = 0
#         bm = b * m % p
#         while S[i][0] != bm:
#             i += 1
#         for j in xrange(i, M, p):
#             l = 1
#             while (S[i][0] + S[i][1] * m) % p ** (l + 1) == 0:
#                 l += 1
#             R[j].append((p, l))

#     # Sieve with algebraic factor base
#     for (r, q) in algebraic:
#         i = 0
#         br = b * r % q
#         while i < len(S) and S[i][0] != br:
#             i += 1
#         for j in xrange(i, M, q):
#             A[j].append((r, q))

#     # Sieve with quadratic factor base
#     return U
