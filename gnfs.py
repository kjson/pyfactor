from sympy import *
import theano 
import math 
import numpy 
from util import *

"""
PS1=$'\e[0;31m$ \e[0m'
"""

def free_parameters(N): 
    """ Returns free parameters used in GNFS """
    # bitlength of N 
    k = math.log(N,2) 
    
    if k >= 110: d = 5.0
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

def algebraic_factor_base(f,bound):
	""" The algebraic for bases in Z[r] where r is a root of f """
	aBase = list()
	assert bound > 2 

	for p in primes_less_than(bound):
		for r in xrange(0,p):
			if f.eval(r) % p == 0:
				aBase.append((r,p))

	return aBase 

def quadratic_character_base(f,l,u):
	""" The quadratic character base of f with primes between l and u """
	qBase = list()

	for p in xrange(l,u):
		if is_probable_prime(p):
			for r in xrange(0,p):
				if f.eval(r) % p == 0:
					qBase.append((r,p))

	return qBase 

def sieve(rBase,aBase,theta,m,b,a):
	""" select m smooth algebraic integers in Z[r] where r is a root of f"""
	pass


def difference_of_squares(smooth_list,N,f):
	""" finds difference of squares """
	pass

def gnfs(N):
	""" Factor given integer N """
    # Select free f, m parameters depending on N 
   	f, m = free_parameters(N)

	# Find rational, algebraic and quadratic factor bases
	rBase = primes_less_than(m)
	aBase = algebraic_factor_base(f,90) # this is brute force 
	qBase = quadratic_character_base(f,90,125) # this is also brute force 

	l = len(rBase) + len(aBase) + len(qBase)
















