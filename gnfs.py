#!/usr/bin/env python

from sympy import *
import theano 
import math 
import numpy  as np
from util import *


# TODO: 
# - Factor bases 
# - Determine N from page 10 of gnfs2.pdf
# - How to determine b from page 9 of gnfs2.pdf
# - Solve U
# - Use gcd to factor N 

""" The General number field sieve """

def free_parameters(N): 
    """ Returns free parameters f,m used in GNFS """

    # bitlength of N 
    k = math.log(N,2) 
    if k >= 110:
    	d = 5.0
    elif 80 < k < 110:
    	d = 4.0
    else:
    	d = 3.0

	# Calc m such that m^d = N 
	m = int(math.ceil(N**(1/d)))

	# Coefficients in base-m rep of N 
	coeffs = base_coeffs(m,N) 

	# Create polynomial 
	x = symbols('x')
	f = poly(x,domain=ZZ) - x

	for i in xrange(0,len(coeffs)):
		f += coeffs[i]*x**i

	return f, m

	# Need speed up with roots mod p
def algebraic_factor_base(f,bound):
	""" Calculates algebraic factor base up to a given bound """
	aBase = list()
	assert bound > 2 

	for p in primes_less_than(bound):
		for r in xrange(0,p):
			if f.eval(r) % p == 0:
				aBase.append((r,p))

	return aBase 

	# Need speed up with roots mod p
def quadratic_character_base(f,l,u):
	""" The quadratic character base of f with primes between l and u """
	qBase = list()

	for p in xrange(l,u):
		if is_probable_prime(p):
			for r in xrange(0,p):
				if f.eval(r) % p == 0:
					qBase.append((r,p))

	return qBase 


def sieve(rBase,aBase,qBase,f,m):
	""" 
	Sieve to find pairs (a,b) such that
		a + b*m is smooth over rBase 
		a + b*theta is smooth over aBase  
		a + b*m is a square in Z and Z[theta]
	""" 

	# Lengths of factor bases 
	lr = len(rBase)
	la = len(aBase)
	lq = len(qBase)
	# Number of columns in matrix 
	lt = lr + la + lq + 1 

	# Number of entries to check
	M = 100

	# Random starting value
	b = random.randrange(1,m)

	# List of M "soil" integers
	S = [(a,b) for a in xrange(0,M)]

	# Lists to keep track of smooth pairs 
	R = [[] for _ in S]
	A = [[] for _ in S]

	# Matrix of all smooth pairs
	U = np.zeros(shape=(lt+1,lt))

	print U

	# Sieve with rational factor base 
	for p in rBase:
		i = 0 
		while (S[i][0] + S[i][1] * m) % p != 0:
			i += 1 

		for j in xrange(i,M,p):
			l = 1 
			while (S[i][0] + S[i][1] * m)% p**(l+1) == 0:
				l +=1 
			R[j].append((p,l))

	# Sieve with algebraic factor base 
	for (r,q) in aBase:
		i = 0
		while S[i][0] + S[i][1]*r % q != 0:
			i+=1

		for j in xrange(i,M,q):
			l=1 
			while S[i][0] + S[i][1] * r % q**(l+1) == 0:
				l+=1
			A[j].append((r,q,l))


	# Check for smoothness in Z and Z[theta]
	for i in xrange(0,len(S)):

		prod = 1
		for (p,l) in R[i]:
			prod *= p**l
		if S[i][0] + S[i][1] * m == prod:
			print " this is Z smooth "

		prod = 1 
		for (r,q,l) in A[i]:
			prod *= q 
		if prod == (-S[i][1])**degree(f) * f.eval(- S[i][0]/S[i][1]):
			print " this is Z[theta] smooth"

	return U


def gnfs(N):
	""" Factor given integer N """

    # Select free f, m parameters depending on N 
    # Should also determine the number of integers to sieve: M
	f, m = free_parameters(N)

	# Find rational, algebraic and quadratic factor bases
	rBase = primes_less_than(m)
	# aBase = algebraic_factor_base(f,m)  
	# qBase = quadratic_character_base(f,m) 

	sieve(rBase,[],[],f,m)























