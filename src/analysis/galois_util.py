# File: galois_util.py
# Author: Christopher Wood

import math

def createBaseElem(base, index, k):
	ii = index
	coeff = []
	for j in range(k):
		coeff.append(ii % base)
		ii = ii / base
	coeff.reverse()
	return coeff

def createExtElem(base, gg, n, m):
	p = []
	n = 2**n
	for j in range (m):
		p.append(createBaseElem(2, gg % n, n))
		gg = gg / n
	p.reverse()
	return p

def coset(s, n):
	order = 2**n - 1
	ns = 0
	for i in range(1, order):
		if (((s * (2 ** i)) % order) == s % order):
			ns = i
			break
	cset = []
	for i in range(ns):
		cset.append((s * (2**i)) % order)
	return cset, min(cset), len(cset)

def cosetLeaders(ss, n):
	leaders = []
	for s in ss:
		cset, leader, size = coset(s, n)
		leaders.append(leader)
	return leaders

def listToString(l):
	s = ""
	for i in l:
		s = s + str(i)
	return s

def toBinString(n, dim):
	bits = []
	for i in range(dim):
		if ((1 << (dim - i - 1)) & n > 0):
			bits.append(1)
		else:
			bits.append(0)
	return bits

def primeFactors(n):
	factors = []
	while (n > 1):
		factor = getSinglePrimeFactor(n)
		if not (factor in factors): 
			factors.append(factor)
		n = n / factor
	return tuple(factors)

def getSinglePrimeFactor(n):
	if (n % 2 == 0):
		return 2
	for d in range(3, int(math.ceil(math.sqrt(n)) + 1), 2):
		if (n % d == 0):
			return d 
	return n

def main():
	print(toBinString(10, 8))
	print(primeFactors(255))
	print(coset(0, 4))
	print(coset(1, 4))
	print(coset(3, 4))
	print(coset(5, 4))
	print(coset(7, 4))
	print(coset(2**4 + 103, 9))

if __name__ == "__main__":
	main()
