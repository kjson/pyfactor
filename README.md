gnfspy
=======

The [General Number Field Sieve](http://en.wikipedia.org/wiki/General_number_field_sieve) in pure Python.

Idea 
-----

There are [many](http://en.wikipedia.org/wiki/General_number_field_sieve#Implementations) implementations of the general number field sieve out there which require excessive parameter choosing to get started. This is designed to be the opposite. For example, 

	>>> gnfs(338954395)
	>>> [5,13,19,274457]

Reading material 
----------------

Most of the math needed is from the following two resources;

- Matthew E. Briggs, [An Introduction to the General Number Field Sieve](http://scholar.lib.vt.edu/theses/public/etd-32298-93111/materials/etd.pdf)

- Carl Pomerance, [A Tale of Two Sieves](http://www.ams.org/notices/199612/pomerance.pdf) 


Todo
-----
- Parallelize with cython

Dependencies
------------

- Python 2.7 or 3.3+

- [Sympy](http://www.sympy.org/)

- [Numpy](http://www.numpy.org/)

- [Cython](http://cython.org/)


