# File: analysis.py 
# Author: Christopher Wood
# Description: Generate all of the metrics for a particular map

# not running as a Sage program, for now.
###  !/usr/bin/env sage -python

import sys
import gc
from math import *
# from guppy import hpy
#from sage.all import *
#from sage.crypto.boolean_function import BooleanFunction
#from sage.crypto.mq.sbox import SBox
from numpy import matrix
import numpy
# sys.path.append("../../GaloisLibrary/")
#from galois import *
#from galois_util import *
from time import *

debug = True

def coset(s, n):
	order = 2**n - 1
	ns = 0
	for i in range(1, order):
		if (((s * (2 ** i)) % order) == s % order):
			ns = i
			break
	cset = []
	for i in range(ns):
		cset.append((s * (2**i)) % order) # this will append s
	return cset, min(cset), len(cset)

# # Simple tests to match the paper
# print(coset(0, 4))
# print(coset(1, 4))
# print(coset(3, 4))
# print(coset(5, 4))
# print(coset(7, 4))

def wt(x):
	bits = len(bin(abs(x))[2:]) # hack O.O
	w = 0
	for i in range(bits):
		if ((x & (1 << i)) > 0):
			w = w + 1
	return w

def testMasks(S):
	order = len(S)
	n = int(log(order, 2))
	for mask in range(1, order): # omit the 0 function
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		bf = BooleanFunction(f)
		print(bf.truth_table(format='int'))

def differentialUniformity(S, aal = -1, aau = -1):
	''' Compute the differential uniformity of the S-box.
	Straight interpretation of the definition.
	'''
	c = 0
	delta = 0
	n = len(S)
	if (aal == -1):
		aal = 1
		aau = n
	for alpha in range(aal, aau): # alpha \not= 0
		for beta in range(n):
			c = 0
			# start = time()
			for z in range(n):
				if ((S[(z ^ alpha)] ^ S[z]) == beta):
					c = c + 1
			if (c > delta): # bump up delta
				delta = c
			# end = time()
			# print(end - start)
	return delta 

def sage_nonlinearity(S, aal = -1, aau = -1):
	''' Nonlinearity of an S-box is the minimum nonlinearity over
	all nonlinear combinations of the coordinate functions.
	'''
	order = len(S)
	n = int(log(order, 2))
	nl = 1 << order # some obnoxious number..
	if (aal == -1):
		aal = 1
		aau = order
	for mask in range(aal, aau): # omit the 0 function
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		bf = BooleanFunction(f)
		bfnl = bf.nonlinearity()
		if (bfnl < nl):
			nl = bfnl # we're after the minimum nonlinearity
	return nl

def nonlinearity(S, aal = -1, aau = -1):
	''' Nonlinearity of an S-box is the minimum nonlinearity over
	all nonlinear combinations of the coordinate functions.
	'''
	global debug
	order = len(S)
	n = int(log(order, 2))
	nl = 1 << order # some obnoxious number..
	if (aal == -1):
		aal = 1
		aau = order
	for mask in range(aal, aau): # omit the 0 function
		if debug:
			print >> sys.stderr, "Mask: " + str(mask)
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		# bf = BooleanFunction(f)
		# bfnl = bf.nonlinearity()
		bfnl = bf_nonlinearity(f, n)
		# print >> sys.stderr, ""
		if (bfnl < nl):
			if debug:
				print >> sys.stderr, "New minimum NL = " + str(bfnl)
			nl = bfnl # we're after the minimum nonlinearity
	return nl

def fwt(f): # fast walsh transform of Boolean function f
	# http://mandala.co.uk/fwt/plugin/Walsh.java
	# http://www.ciphersbyritter.com/ARTS/MEASNONL.HTM
	wf = []
	for x in f:
		if x == 0:
			wf.append(1)
		elif x == 1:
			wf.append(-1)
		else:
			raise Exception("Invalid truth table.")
	k = len(f) # k = 2^n
	n = int(log(k, 2))
	sw = 1
	bs = k - 1
	while True:
		li = 0
		bs = bs >> 1
		for b in range(bs, -1, -1):
			ri = li + sw
			for p in range(0, sw):
				a = wf[li]
				b = wf[ri]
				wf[li] = a + b
				wf[ri] = a - b
				li = li + 1
				ri = ri + 1
			li = ri
		sw = (sw << 1) & (k - 1)
		if (sw == 0):
			break
	return wf

def bf_nonlinearity(f, n):
	fw = fwt(f)
	for i in range(len(fw)):
		fw[i] = abs(fw[i])
	return ((2**(n-1)) - (max(fw) / 2))

def sage_differentialApproximationTable(S, formatOutput = False):
	s = mq.SBox(S)
	order = len(S)
	diffMatrix = s.difference_distribution_matrix()
	maxDiff = 0
	for i in range(order):
		for j in range(order):
			if (diffMatrix[i][j] > maxDiff):
				maxDiff = diffMatrix[i][j]

	if (formatOutput):
		result = ""
		for i in range(order):
			for j in range(order):
				result = result + str(diffMatrix[i][j]) + ","
		result = result + "\n" + str(maxDiff) # output the maximum difference too...
		return result
	else:
		return diffMatrix, maxDiff

def differenceDistributionTable(S, formatOutput = False):
	dimension = len(S)
	ndd = numpy.zeros((dimension, dimension), dtype=numpy.uint8)
	maxDiff = 0
	for x1 in range(dimension):
		for x2 in range(dimension):
			dx = x1 ^ x2
			dy = S[x1] ^ S[x2]
			ndd[(dx,dy)] = ndd[(dx,dy)] + 1
			if ndd[(dx,dy)] > maxDiff:
				maxDiff = ndd[(dx,dy)]

	if (formatOutput):
		result = ""
		for i in range(dimension):
			for j in range(dimension):
				result = result + str(ndd[i][j]) + ","
		result = result + "\n" + str(maxDiff) # output the maximum difference too...
		return result
	else:
		return ndd, maxDiff

def branchNumber(S):
	order = len(S)
	n = int(log(order, 2))
	bn = 1 << order # some obnoxious number..
	for a in range(0, order):
		for b in range(0, order):
			if (a != b):
				twt = wt(a ^ b) + wt(S[a] ^ S[b])
				if (twt < bn):
					bn = twt
	return bn

def linearApproximationTable(S, formatOutput = False):
	NL = {}
	bias = {}
	dimension = len(S)
	n = int(log(dimension, 2))
	for i in range(dimension):
		for j in range(dimension):
			NL[(i,j)] = 0
	computeNL(S, int(n), len(S), NL)

	# Put NL in matrix form
	nlm = numpy.zeros((dimension, dimension), dtype=numpy.uint8)
	maxEntry = 0
	for entry in NL.keys():
		nlm[entry[0]][entry[1]] = NL[entry] - (dimension / 2)
		if (abs(nlm[entry[0]][entry[1]]) > maxEntry):
			maxEntry = abs(nlm[entry[0]][entry[1]])

	if (formatOutput):
		result = ""
		for i in range(dimension):
			for j in range(dimension):
				result = result + str(NL[(i, j)]) + ","
		result = result + "\n" + str(maxEntry)
		return result
	else:
		return nlm, maxEntry # maxEntry is maximum absolute bias

def isBalanced(S): 
	''' Balanced == permutation. Simple check.
	'''
	order = len(S)
	n = int(log(order, 2))
	seen = []
	for x in range(order):
		if x in seen:
			return False
		seen.append(x)
	return seen

def sac(S, formatOutput = False):
	order = len(S)
	n = int(log(order, 2))
	bitBucket = []
	for i in range(n):
		#bitBucket.append()
		bucket = []
		for j in range(n):
			bucket.append(0)
		bit = 1 << i 
		for x in range(order):
			xorDiff = S[x] ^ S[x ^ bit]
			for b in range(n):
				if (((1 << b) & xorDiff) != 0):
					bucket[b] = bucket[b] + 1	
			#bitBucket[i] = bitBucket[i] + ((S[x] ^ S[x ^ bit]) % 2)
		bitBucket.append(bucket)
	if (formatOutput):
		result = ""
		for b in bitBucket:
			result = result + str(b) + ","
		return result
	else:
		return bitBucket

def algebraicImmunityComponent(S): 
	''' TODO : from Claret chapter on vectorial boolean functions (page 14)
	'''
	order = len(S)
	n = int(log(order, 2))
	opt_bfr = 1 << n # some large number...
	for mask in range(1, order): # omit the 0 function
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		bf = BooleanFunction(f)
		bfr = bf.algebraic_immunity()
		if (bfr < opt_bfr):
			opt_bfr = bfr
		del bf
		gc.collect()
	return opt_bfr

def algebraicImmunityBiAffine(s, a):
	''' Gamma = (t/s)^ceil(t/r), t = s(s + 2) + 1
		Gamma' = ((t-r)/s)^ceil((t-r)/s)
		r = # biaffine equations for s/a
	'''
	t = s * (s + 2) + 1
	r = 0
	if (a == (2**s - 2)):
		r = (3*s) - 1
	else:
		r = numLinearIndependentBiAffineEqns(a, s)
	if r == 0:
		return 0, 0
	t = float(t)
	r = float(r)
	s = float(s)
	gamma = (t / s)**(ceil(t / r))
	gammap = ((t - r) / s)**(ceil((t - r) / s)) 
	return gamma, gammap

def algebraicImmunityQuadratic(s, a):
	''' Gamma = (t/s)^ceil(t/r), t = s(2s + 1) + 1
		Gamma' = ((t-r)/s)^ceil((t-r)/s)
		r = # quadratic equations for s/a
	'''
	t = s * (2 * s + 1) + 1
	r = 0
	if (a == (2**s - 2)):
		r = (5*s) - 1
	else:
		r = numLinearIndependentQuadraticEqns(a, s)
	if r == 0:
		return 0, 0
	t = float(t)
	r = float(r)
	s = float(s)
	gamma = (t / s)**(ceil(t / r))
	gammap = ((t - r) / s)**(ceil((t - r) / s)) 
	return gamma, gammap

def numLinearIndependentQuadraticEqns(a, n):
	if (n % 2 != 0):
		raise Exception("Invalid parameters n/m")
	m = n / 2

	# Step 1
	cst_l = [] # coset leaders 
	cst_s = [] # coset sizes

	# a coset leader/size
	cset, leader, size = coset(a, n)
	cst_l.append(leader)
	cst_s.append(size)

	# First value...
	for k in range(n):
		ak = (2**k + a) # % (2**n) # ADDED MOD
		cset, leader, size = coset(ak, n)
		cst_l.append(leader)
		cst_s.append(size)

	# Second value...
	for k in range(1, m):
		ak = ((2**k + 1) * a) # % (2**n) # ADDED MOD
		cset, leader, size = coset(ak, n)
		cst_l.append(leader)
		cst_s.append(size)

	# Step 2: compute the special value
	ak = ((2**m + 1)*a) # % (2**n) # ADDED MOD
	cset, leader, size = coset(ak, n) # CAW: n or m? Not specified in the algorithm.
	cstm_l = leader
	cstm_s = size

	# Total number of elements = n + (m - 1) + 1 = n + m
	N = n + m

	# Step 3: perform sort
	swapped = True
	while swapped:
		swapped = False
		for i in range(N - 1):
			if (cst_l[i] > cst_l[i + 1]):
				j = i + 1
				tmpl = cst_l[i]
				cst_l[i] = cst_l[j]
				cst_l[j] = tmpl
				tmps = cst_s[i]
				cst_s[i] = cst_s[j]
				cst_s[j] = tmps
				swapped = True

	# Step 4: Initialize k and eqn (# of eqns)
	k = 0
	eqn = 0

	# Step 5: Iterate (longer paper version)
	while (k < N): 
		cstlk = wt(cst_l[k])
		if (0 < cstlk and cstlk <= 2):
			eqn = eqn + n
		else:
			if (cst_s[k] < n):
				eqn = eqn + (n - cst_s[k])
			if (k != N - 1 and cst_l[k] == cst_l[k + 1]):
				eqn = eqn + cst_s[k]
		k = k + 1

	# Line 11
	cstlmlwt = wt(cstm_l)
	if (0 < cstlmlwt and cstlmlwt <= 2):
		eqn = eqn + m
	else:
		if (cstm_s < m):
			eqn = eqn + (m - cstm_s)
		k = 0
		while (k < N):
			if (cstm_l == cst_l[k]):
				# eqn = eqn + cstm_s
				break
			k = k + 1

	# Finally, return the equation count...
	return eqn

def numLinearIndependentBiAffineEqns(a, n):
	''' Number of biaffine equations for y = x^a over F_2^n
	'''
	cst_l = [] # coset leaders 
	cst_s = [] # coset sizes

	# a coset leader/size
	cset, leader, size = coset(a, n)
	cst_l.append(leader)
	cst_s.append(size)

	for k in range(n):
		ak = ((2**k) + a) # % (2**n) # ADDED MOD
		cset, leader, size = coset(ak, n)
		cst_l.append(leader)
		cst_s.append(size)

	# Total number of items = n + 1
	N = n + 1

	# Perform sort
	swapped = True
	while swapped:
		swapped = False
		for i in range(N - 1):
			if (cst_l[i] > cst_l[i + 1]):
				j = i + 1
				tmpl = cst_l[i]
				cst_l[i] = cst_l[j]
				cst_l[j] = tmpl
				tmps = cst_s[i]
				cst_s[i] = cst_s[j]
				cst_s[j] = tmps
				swapped = True

	# Now actually do the counting...
	eqn = 0
	k = 0
	while (k < N):
		if (wt(cst_l[k]) == 0):
			eqn = eqn + (n - 1)
		elif (wt(cst_l[k]) == 1):
			eqn = eqn + n
		else:
			if (cst_s[k] < n):
				eqn = eqn + (n - cst_s[k])
			if (k != N - 1 and cst_l[k] == cst_l[k+1]):
				eqn = eqn + cst_s[k]
		k = k + 1

	return eqn

def correlationImmunity(S):
	order = len(S)
	n = int(log(order, 2))
	for t in range(n):
		if isCorrelationImmune(S, order, n, t) == False:
			return t - 1
	raise Exception("ERROR: CI SHOULD BE IN [0,N]")

def sage_isCorrelationImmune(S, order, n, t):
	for mask in range(1, order): # omit the 0 function
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		bf = BooleanFunction(f)
		bfci = bf.correlation_immunity()
		print >> sys.stderr, bfci
		if (bfci < t):
			return False
		del(bf)
	return True

def isCorrelationImmune(S, order, n, t):
	for mask in range(1, order): # omit the 0 function
		f = []
		for x in range(0, order):
			s = 0
			for i in range(0, n):
				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
					s = s ^ 1
			f.append(s)
		spectrum = fwt(f)
		for x in range(order):
			if (wt(x) <= t and spectrum[x] != 0):
				return False
	return True

def resiliency(S):
	order = len(S)
	n = int(log(order, 2))
	for t in range(n):
		for mask in range(1, order): # omit the 0 function
			f = []
			for x in range(0, order):
				s = 0
				for i in range(0, n):
					if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
						s = s ^ 1
				f.append(s)
			spectrum = fwt(f)
			if (spectrum[0] == 0): # balanced
				for u in range(order):
					if (wt(u) <= t and spectrum[u] != 0):
						return t - 1
			else:
				return t - 1

# def resiliency(S):
# 	order = len(S)
# 	n = int(log(order, 2))
# 	for t in range(n):
# 		if isResilient(S, order, n, t) == False:
# 			return t - 1
# 	raise Exception("ERROR: RESILIENCY SHOULD BE IN [0,N]")

# def isResilient(S, order, n, t):
# 	for mask in range(1, order): # omit the 0 function
# 		f = []
# 		for x in range(0, order):
# 			s = 0
# 			for i in range(0, n):
# 				if ((mask & (1 << i)) > 0) and ((S[x] & (1 << i)) > 0):
# 					s = s ^ 1
# 			f.append(s)
# 		bf = BooleanFunction(f)
# 		bfr = bf.resiliency_order()
# 		if (bfr < t):
# 			return False
# 		del(bf)
# 	return True

# Compute and return the input/output bitwise-addition sum.
def ioSum(a, b, x, y, bits):
	newX = a & x
	newY = b & y
	xw = wt(newX) % 2 # XOR of bits
	yw = wt(newY) % 2 # XOR of bits
	return (xw ^ yw) % 2 # mod 2 isn't really necessary

# Compute the N_L (linear approximation) table 
def computeNL(S, bits, dimension, NL):
	count = 0
	for a in range(dimension):
		for b in range(dimension):
			for x in range(dimension):
				for y in range(dimension):
					if ((S[x] == y) and (ioSum(a, b, x, y, bits) % 2 == 0)): # mod 2
						if not ((a,b) in NL.keys()):
							NL[(a,b)] = 0
						NL[(a,b)] = NL[(a,b)] + 1 # bump it up
			count = count + 1

def computeBias(NL, bias, dimension):
	for a in range(dimension):
		for b in range(dimension):
			if not ((a,b) in bias.keys()):
				bias[(a,b)] = 0
			bias[(a,b)] = (float(NL[(a,b)]) - (float(dimension) / 2)) / dimension

def components(x):
    return [i for i in range(x + 1) if x & i == i]

def genANF(sbox):
    num = len(sbox)
    dim = int(log(num, 2))
    if 2 ** dim != num:
        raise InputError # make sure we're a power of 2 (clearly...)

    ls = []
    for d in range(dim):
        dc = []
        for n in range(num):
            sign = 0
            for i in components(n): # the support of the BF
                sign = sign ^ ((sbox[i] >> d) & 0x1)
            dc.append(sign)
        ls.append(dc)
    return ls

def genTT(sbox):
	num = len(sbox)
	dim = int(log(num, 2))
	if 2 ** dim != num:
		raise InputError # make sure we're a power of 2 (clearly...)
	tt = []
	for d in range(dim):
		tts = []
		for i in range(num):
			if ((sbox[i] >> d) & 0x1):
				tts.append(1)
			else:
				tts.append(0)
		tt.append(tts)
	return tt

def test():
	global debug
	
	# S = [0x1, 0x2, 0x3, 0x0]
	# testMasks(S) # A simple test...
	# print(sac(S))
	# S = [0xE, 0x4, 0xD, 0x1, 0x2, 0xF, 0xB, 0x8, 0x3, 0xA, 0x6, 0xC, 0x5, 0x9, 0x0, 0x7]
	# S = [7,6,0,4,2,5,1,3]
	# print(sage_differentialApproximationTable(S, True))
	# print(differenceDistributionTable(S, True))
	# print(linearApproximationTable(S, True))
	# test
	# print(numLinearIndependentBiAffineEqns(254, 8))
	# print(numLinearIndependentQuadraticEqns(254, 8))

	print >> sys.stderr, "Beginning analysis tests..."

	# Power mapings with zero biaffine equations
	print >> sys.stderr, "Testing: Power mapings with zero biaffine equations"
	assert numLinearIndependentBiAffineEqns(11, 8)  == 0
	assert numLinearIndependentBiAffineEqns(23, 8)  == 0
	assert numLinearIndependentBiAffineEqns(29, 8)  == 0
	assert numLinearIndependentBiAffineEqns(61, 8)  == 0
	assert numLinearIndependentBiAffineEqns(13, 10) == 0
	assert numLinearIndependentBiAffineEqns(19, 10) == 0
	assert numLinearIndependentBiAffineEqns(23, 10) == 0
	assert numLinearIndependentBiAffineEqns(43, 10) == 0
	assert numLinearIndependentBiAffineEqns(47, 10) == 0
	assert numLinearIndependentBiAffineEqns(53, 10) == 0
	assert numLinearIndependentBiAffineEqns(59, 10) == 0
	assert numLinearIndependentBiAffineEqns(61, 10) == 0
	assert numLinearIndependentBiAffineEqns(71, 10) == 0
	assert numLinearIndependentBiAffineEqns(79, 10) == 0
	assert numLinearIndependentBiAffineEqns(89, 10) == 0
	assert numLinearIndependentBiAffineEqns(109, 10) == 0
	assert numLinearIndependentBiAffineEqns(119, 10) == 0
	assert numLinearIndependentBiAffineEqns(125, 10) == 0
	assert numLinearIndependentBiAffineEqns(151, 10) == 0
	assert numLinearIndependentBiAffineEqns(175, 10) == 0
	assert numLinearIndependentBiAffineEqns(191, 10) == 0
	assert numLinearIndependentBiAffineEqns(221, 10) == 0
	assert numLinearIndependentBiAffineEqns(245, 10) == 0
	assert numLinearIndependentBiAffineEqns(251, 10) == 0
	print >> sys.stderr, "Passed."

	# Power mapings with zero quadratic equations
	print >> sys.stderr, "Testing: Power mapings with zero quadratic equations"
	assert numLinearIndependentQuadraticEqns(27, 8) == 0
	assert numLinearIndependentQuadraticEqns(27, 10) == 0
	assert numLinearIndependentQuadraticEqns(45, 10) == 0
	assert numLinearIndependentQuadraticEqns(51, 10) == 0
	assert numLinearIndependentQuadraticEqns(53, 10) == 0
	assert numLinearIndependentQuadraticEqns(75, 10) == 0
	assert numLinearIndependentQuadraticEqns(87, 10) == 0
	assert numLinearIndependentQuadraticEqns(105, 10) == 0
	assert numLinearIndependentQuadraticEqns(111, 10) == 0
	assert numLinearIndependentQuadraticEqns(117, 10) == 0
	assert numLinearIndependentQuadraticEqns(123, 10) == 0
	assert numLinearIndependentQuadraticEqns(183, 10) == 0
	assert numLinearIndependentQuadraticEqns(237, 10) == 0
	assert numLinearIndependentQuadraticEqns(251, 10) == 0
	print >> sys.stderr, "Passed."

	print >> sys.stderr, "Testing: Highly nonlinear power mapping tests"
	assert numLinearIndependentBiAffineEqns(31, 8)    == 16
	assert numLinearIndependentQuadraticEqns(31, 8)   == 36
	assert numLinearIndependentBiAffineEqns(91, 8)    == 16
	assert numLinearIndependentQuadraticEqns(91, 8)   == 36
	assert numLinearIndependentBiAffineEqns(127, 8)   == 23
	assert numLinearIndependentQuadraticEqns(127, 8)  == 39
	assert numLinearIndependentBiAffineEqns(5, 10)   == 10
	assert numLinearIndependentQuadraticEqns(5, 10)  == 40
	assert numLinearIndependentBiAffineEqns(13, 10)   == 0
	assert numLinearIndependentQuadraticEqns(13, 10)  == 20
	assert numLinearIndependentBiAffineEqns(17, 10)   == 15
	assert numLinearIndependentQuadraticEqns(17, 10)  == 40
	assert numLinearIndependentBiAffineEqns(25, 10)   == 5
	assert numLinearIndependentQuadraticEqns(25, 10)  == 10
	assert numLinearIndependentBiAffineEqns(41, 10)   == 5
	assert numLinearIndependentQuadraticEqns(41, 10)  == 5
	assert numLinearIndependentBiAffineEqns(49, 10)   == 5
	assert numLinearIndependentQuadraticEqns(49, 10)  == 15
	assert numLinearIndependentBiAffineEqns(79, 10)   == 0
	assert numLinearIndependentQuadraticEqns(79, 10)  == 20
	assert numLinearIndependentBiAffineEqns(107, 10)   == 5
	assert numLinearIndependentQuadraticEqns(107, 10)  == 15
	assert numLinearIndependentBiAffineEqns(181, 10)  == 15
	assert numLinearIndependentQuadraticEqns(181, 10) == 35
	assert numLinearIndependentBiAffineEqns(205, 10)   == 10
	assert numLinearIndependentQuadraticEqns(205, 10)  == 40
	assert numLinearIndependentBiAffineEqns(511, 10)   == 29
	assert numLinearIndependentQuadraticEqns(511, 10)  == 49
	print >> sys.stderr, "Passed."

	print >> sys.stderr, "Testing: Cryptographically significant power mappings"
	assert numLinearIndependentBiAffineEqns(29, 8)    == 0
	assert numLinearIndependentQuadraticEqns(29, 8)   == 24
	assert numLinearIndependentBiAffineEqns(79, 10)    == 0
	assert numLinearIndependentQuadraticEqns(79, 10)   == 20
	assert numLinearIndependentBiAffineEqns(223, 10)    == 5
	assert numLinearIndependentQuadraticEqns(223, 10)   == 5	
	assert numLinearIndependentBiAffineEqns(731, 12)    == 9
	assert numLinearIndependentQuadraticEqns(731, 12)   == 9
	assert numLinearIndependentBiAffineEqns(319, 14)    == 0
	assert numLinearIndependentQuadraticEqns(319, 14)   == 28
	assert numLinearIndependentBiAffineEqns(1883, 14)    == 0
	assert numLinearIndependentQuadraticEqns(1883, 14)   == 0
	assert numLinearIndependentBiAffineEqns(73, 12)    == 9
	assert numLinearIndependentQuadraticEqns(73, 12)   == 9
	assert numLinearIndependentBiAffineEqns(341, 12)    == 10
	assert numLinearIndependentQuadraticEqns(341, 12)   == 10
	assert numLinearIndependentBiAffineEqns(731, 12)    == 9
	assert numLinearIndependentQuadraticEqns(731, 12)   == 9
	assert numLinearIndependentBiAffineEqns(853, 12)    == 10
	assert numLinearIndependentQuadraticEqns(853, 12)   == 10
	print >> sys.stderr, "Passed."

	print >> sys.stderr, "Select test cases"
	assert numLinearIndependentBiAffineEqns(127, 8)    == 23
	assert numLinearIndependentQuadraticEqns(127, 8)   == 39
	assert numLinearIndependentBiAffineEqns(191, 8)    == 23
	assert numLinearIndependentQuadraticEqns(191, 8)   == 39
	assert numLinearIndependentBiAffineEqns(223, 8)    == 23
	assert numLinearIndependentQuadraticEqns(223, 8)   == 39
	assert numLinearIndependentBiAffineEqns(239, 8)    == 23
	assert numLinearIndependentQuadraticEqns(239, 8)   == 39
	assert numLinearIndependentBiAffineEqns(247, 8)    == 23
	assert numLinearIndependentQuadraticEqns(247, 8)   == 39
	assert numLinearIndependentBiAffineEqns(251, 8)    == 23
	assert numLinearIndependentQuadraticEqns(251, 8)   == 39
	assert numLinearIndependentBiAffineEqns(253, 8)    == 23
	assert numLinearIndependentQuadraticEqns(253, 8)   == 39
	assert numLinearIndependentBiAffineEqns(254, 8)    == 23
	assert numLinearIndependentQuadraticEqns(254, 8)   == 39
	print >> sys.stderr, "Passed."

	# # HERE NOW MKAY
	print >> sys.stderr, "Testing: FWT"
	assert fwt([0,1,1,1,0,1,0,0]) == [0, 4, 0, 4, -4, 0, 4, 0]
	n = 3
	f = [1,0,0,1,1,1,0,0]
	# print(f)
	assert fwt(f) == [0, 0, -4, -4, 0, 0, 4, -4]
	print >> sys.stderr, "Passed."

	print >> sys.stderr, "Testing: Nonlinearity calculation"
	assert bf_nonlinearity(f, n) == 2
	print >> sys.stderr, "Passed."

def computeAll(S, n, a):

	print >> sys.stderr, "Computing all for n = " + str(n) + ", a = " + str(a)
	print >> sys.stderr, S

	vals = []
	# v = differentialUniformity(S)
	# vals.append(v)
	# print >> sys.stderr, v
	print >> sys.stderr, "Computing nonlinearity..."
	v = nonlinearity(S)
	vals.append(v)
	print >> sys.stderr, v
	print >> sys.stderr, "Computing branch number..."
	v = branchNumber(S)
	vals.append(v)
	print >> sys.stderr, v
	print >> sys.stderr, "Computing algebraic immunity..."
	val = algebraicImmunityComponent(S) # compute AI in R, not using SAGE
	vals.append(v)
	print >> sys.stderr, val
	print >> sys.stderr, "Computing alg immunity biaffine..."
	v = algebraicImmunityBiAffine(n, a)
	vals.append(v)
	print >> sys.stderr, v
	print >> sys.stderr, "Computing alg immunity quadratic..."
	v = algebraicImmunityQuadratic(n, a)
	vals.append(v)
	print >> sys.stderr, v
	print >> sys.stderr, "Computing number biaffine..."
	v = numLinearIndependentBiAffineEqns(a, n)
	vals.append(v)
	print >> sys.stderr, v
	print >> sys.stderr, "Computing number quadratic..."
	v = numLinearIndependentQuadraticEqns(a, n)
	vals.append(v)
	print >> sys.stderr, v
	v = correlationImmunity(S)
	vals.append(v)
	print >> sys.stderr, v
	v = resiliency(S)
	vals.append(v)
	print >> sys.stderr, v
	v = sac(S, formatOutput = True)
	vals.append(v)
	print >> sys.stderr, v
	# val = differenceDistributionTable(S, formatOutput = True)
	# print >> sys.stderr, val
	# vals.append(val)
	# val = linearApproximationTable(S, formatOutput = True)
	# print >> sys.stderr, val
	# vals.append(val)

	return vals

def parseSinglePoly(n, p):
	if p == "0": # special case
		return 0, "0x0"
	data = p.split(" + ")
	bs = []
	for i in range(n):
		bs.append(0)
	for b in data:
		if (b == '1'):
			bs[n - 1] = 1
		elif (b == 'x'):
			bs[n - 2] = 1
		else:
			bs[n - int(b[2:]) - 1] = 1	
	bsi = ""
	for b in bs:
		bsi = bsi + str(b)
	val = int(bsi, 2)
	return val, str(hex(int(bsi, 2)))

def parsePolynomialFile(fname):
	f = open(fname, 'r')
	polynomial = f.readline().strip()
	header = f.readline().strip().split(" ")
	n, d = int(header[0]), int(header[1])

	# Read in the map, do the conversion, etc etc
	domain = []
	e = f.readline().strip()
	while (len(e) > 0):
		ie = e
		e = f.readline().strip()
		oe = e
		iei, ieh = parseSinglePoly(n, ie)
		oei, oeh = parseSinglePoly(n, oe)
		domain.append((iei, oeh))

		# Advance...
		e = f.readline().strip()

	# Sort the domain map by the input values (their integer indices)
	swapped = True
	while swapped:
		swapped = False
		for i in range(len(domain)):
			for j in range(i + 1, len(domain)):
				if (i != j):
					if (domain[i][0] > domain[j][0]):
						tmp = domain[i]
						domain[i] = domain[j]
						domain[j] = tmp
						swapped = True

	sortedDomain = []
	for i in range(len(domain)):
		sortedDomain.append(int(domain[i][1], 0))
	return sortedDomain, n, d, polynomial

def parseHexFile(fname):
	f = open(fname, 'r')
	polynomial = f.readline().strip() # strip out the polynomial
	header = f.readline().strip() # strip out the actual header
	sa = header.split(" ")
	n = int(sa[0].strip()) # order of field
	a = int(sa[1].strip()) # exponent for power mapping
	chunk = f.readline()
	data = chunk.split(",")
	elems = []
	for i in range(len(data)):
		if ('x' in data[i]):
			elems.append(int(data[i], 0))
	return elems, n, a, polynomial

def main():
	if len(sys.argv) == 1: # use this one if no cmd line arguments given
		test()
	elif len(sys.argv) == 2:
		fnames = open(sys.argv[1], 'r')
		for fn in fnames:
			S, n, a, polynomial = parseHexFile(fn.strip())
			start = time()
			vals = computeAll(S, n, a)
			end = time()

			print >> sys.stderr, "Elapsed time: " + str(end - start)
			print(polynomial)
			print(str(n) + " " + str(a))
			for v in vals:
				print(v)
	elif len(sys.argv) == 5:
		S, n, a, polynomial = parseHexFile(sys.argv[1])
		mode = int(sys.argv[2])
		if mode == 0:
			start = time()
			vals = computeAll(S, n, a)
			end = time()

			print >> sys.stderr, "Elapsed time: " + str(end - start)
			print(polynomial)
			print(str(n) + " " + str(a))
			for v in vals:
				print(v)
		elif mode == 1: # DUM
			val = differentialUniformity(S, aal = int(sys.argv[3]), aau = int(sys.argv[4]))
			print(val)
		elif mode == 2: # NL: python analysis.py p n 2 aal aau
			val = nonlinearity(S, aal = int(sys.argv[3]), aau = int(sys.argv[4]))
			print(val)
		elif mode == 3: # BN
			val = 0
			print(val)
		elif mode == 4: # AI_BIAFFINE
			val = 0
			print(val)
		elif mode == 5: # AI_QUAD
			val = 0
			print(val)
		elif mode == 6: # CORRELATION_IMMUNITY
			val = 0
			print(val)
		elif mode == 7: # RESILIENCY
			val = 0
			print(val)
		elif mode == 8: # SAC
			val = 0
			print(val)
	elif len(sys.argv) == 8:
		# NL: python analysis.py p n ip d 2 aal aau
		# 16-bit S-box: python analysis.py 2 16 10000000000101101 65534 2 1 65526
		p = int(sys.argv[1])
		n = int(sys.argv[2])
		bs = sys.argv[3]
		d = int(sys.argv[4])
		mode = int(sys.argv[5])
		lb = int(sys.argv[6])
		ub = int(sys.argv[7])
		h = hpy()

		# Create the field and the mapping given the particular exponent d (this may take time)
		ip = GFElem([], binString = bs)
		field = GF(p, n, ip)
		S = []
		for e in field.getElems():
			S.append(field.power(e, d).toInt(n))

		if mode == 1:
			val = differentialUniformity(S, aal = lb, aau = ub)
			print(val)
		elif mode == 2:
			val = nonlinearity(S, aal = lb, aau = ub)
			print(val)
		print h.heap()
	else:
		print "Usage: %s <bf_file>" % sys.argv[0]
		print "Analyze the properties of the input (n,n) Boolean function"
		sys.exit(1)

# Entry point
if __name__ == "__main__":
	main()
