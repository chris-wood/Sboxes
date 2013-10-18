# File: galois.py
# Author: Christopher Wood
# NOTE: all polynomials are stored as little endian lists of coefficients

# For object serialization
import pickle

# For profiling purposes
import time
import sys

# Exception handling
import traceback

# For external, friendly utilities that are "galois-free" - not coupled to these classes 
from galois_util import *

class GFElem:
	def __init__(self, coefficients, binString = None, id = -1, base = -1, exp = -1): 
		if (id != -1): # override coefficients
			p = []
			gg = id
			for j in range(exp):
				p.append(gg % base)
				gg = gg / base
			p.reverse()
			coefficients = p
		if (binString != None):
			coefficients = []
			for i in range(len(binString)):
				coefficients.append(int(binString[i]))

		# Create the coefficient collection, pruning leading zero degrees if necessary
		self.coeff = []
		coeff = self.prune(coefficients[:])
		for i in range(len(coeff) - 1, -1, -1):
			self.coeff.append(coeff[i])

	def toBinary(self, n): # assumes characteristic 2
		bits = []
		degree = len(self.coeff)
		for i in range(n):
			if (i < degree):
				bits.append(self.coeff[i])
			else:
				bits.append(0)
		bits.reverse()
		return bits

	def toInt(self, n): # assumes characteristic 2
		bits = self.toBinary(n)
		val = 0
		for i in range(len(bits)):
			if (bits[i] == 1):
				val = val | (1 << (n - i - 1))
		return val

	def prune(self, coefficients):
		foundNonZero = False
		index = 0
		while (not foundNonZero and index < len(coefficients)): 
			if (coefficients[index] != 0): 
				foundNonZero = True
			else:
				coefficients.pop(0)
		return coefficients

	def degree(self):
		d = -1
		for i in range(len(self.coeff)):
			if (self.coeff[i] != 0):
				d = i
		return d

	def scalarMult(self, s):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(self.coeff[i] * s)
		coeff.reverse()
		return GFElem(coeff)

	def shiftLeft(self, bits):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(self.coeff[i])
		for i in range(bits):
			coeff.append(0)
		for i in range(bits):
			for j in range(len(coeff) - 1, 0, -1):
				coeff[j] = coeff[j - 1] # shifting to the left -> we have little endian order!
			coeff[0] = 0
		coeff.reverse()
		return GFElem(coeff)

	def shiftRight(self, bits):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(self.coeff[i])
		for i in range(bits):
			for j in range(len(coeff)):
				if (j < len(coeff) - 1):
					coeff[j] = coeff[j + 1] # shifting to the left -> we have little endian order!
					coeff[j + 1] = 0
		coeff.reverse()
		return GFElem(coeff)

	def isUnit(self):
		if (len(self.coeff) == 0):
			return False
		elif (self.coeff[0] == 1 and self.degree() == 0):
			return True
		else:
			return False

	def isZero(self):
		for i in range(len(self.coeff)):
			if (self.coeff[i] != 0):
				return False
		return True

	def getCoeff(self, n = -1):
		c = self.coeff[:]
		if (n != -1):
			if (len(c) < n):
				for i in range(n - len(self.coeff)):
					c.append(0)
		return c

	def bin(self, fill = False, n = 0):
		result = ""
		for c in self.coeff:
			result = result + str(c)
		if (fill == True):
			result = result + ("0" * (n - len(result)))
		return result[::-1]

	def copy(self):
		c = self.coeff[:]
		c.reverse()
		return GFElem(c)

	def lsb(self, i):
		if (i < len(self.coeff)):
			return self.coeff[i]
		else:
			raise TypeError() # big big problem.

	def wt(self):
		count = 0
		for c in self.coeff:
			if (c > 0):
				count = count + c
		return count

	def nterms(self):
		count = 0
		for c in self.coeff:
			if (c > 0):
				count = count + 1
		return count

	def __eq__(self, other):
		if isinstance(other, GFElem):
			match = True
			length = max(len(self.coeff), len(other))
			for i in range(length):
				if (self[i] != other[i]):
					match = False
					break
			return match

	def __len__(self):
		return len(self.coeff)

	def __getitem__(self, key):
		if (len(self.coeff) <= key):
			for i in range(key - len(self.coeff) + 1):
				self.coeff.append(0)
		return self.coeff[key]

	def __setitem__(self, key, value):
		if (len(self.coeff) <= key):
			for i in range(key - len(self.coeff) + 1):
				self.coeff.append(0)
		self.coeff[key] = value

	def toString(self, trim = False):
		if (trim == True):
			poly = ""
			self.coeff.reverse()
			exp = len(self.coeff) - 1
			for x in range(len(self.coeff)):
				if (self.coeff[x] > 0):
					if (str(exp) == "0"):
						poly = poly + "1 + "	
					elif (str(exp) == "1"):
						poly = poly + "x + "
					else:
						poly = poly + "x^" + str(exp) + " + "
				exp = exp - 1
			poly = poly[0:len(poly) - 3] + ""
			if (poly == ""):
				poly = "0"
			self.coeff.reverse()
			return poly
		else:
			return str(self)

	def __str__(self):
		poly = "("
		self.coeff.reverse()
		exp = len(self.coeff) - 1
		for x in range(len(self.coeff)):
			if (self.coeff[x] > 0):
				poly = poly + str(self.coeff[x]) + "x^" + str(exp) + " + "
			exp = exp - 1
		poly = poly[0:len(poly) - 3] + ")"
		if (poly == ")"):
			poly = "0"
		self.coeff.reverse()
		return poly

class GFExtensionElem:
	def __init__(self, coefficients, binString = None, smallSize = 0):
		if (binString != None and smallSize > 0):
			coefficients = []
			if (len(binString) % 2 != 0):
				raise Exception("Invalid binary string in GFExtensionElem constructor")
			num = len(binString) / smallSize
			for i in range(num):
				coefficients.append(GFElem([], binString = binString[smallSize*i:smallSize*i + smallSize]))
		self.coeff = []
		coefficients = self.prune(coefficients[:])
		for i in range(len(coefficients) - 1, -1, -1):
			self.coeff.append(coefficients[i])

	def toBinary(self, n, m): # GF((2^n)^m)
		bits = []
		coeff = []
		for c in self.coeff:
			coeff.append(c)
		coeff.reverse()
		bitList = []
		for i in range(m):
			if (i < len(self.coeff)):
				subBits = self.coeff[i].toBinary(n)
				bitList.append(subBits)
			else:
				fake = GFElem([])
				subBits = fake.toBinary(n)
				bitList.append(subBits)
		bitList.reverse()
		for l in bitList:
			for b in l:
				bits.append(b)
		return bits

	def prune(self, coefficients):
		foundNonZero = False
		index = 0
		zero = GFElem([])
		while (not foundNonZero and index < len(coefficients)): 
			if (not coefficients[index].isZero()): 
				foundNonZero = True
			else:
				coefficients.pop(0)
		return coefficients

	def degree(self):
		d = -1
		for i in range(len(self.coeff)):
			if (not self.coeff[i].isZero()):
				d = i
		return d

	def scalarMult(self, s, field):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(field.g_mult(self.coeff[i], s))
		coeff.reverse()
		return GFExtensionElem(coeff)

	def shiftLeft(self, bits):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(self.coeff[i])
		for i in range(bits):
			coeff.append(GFElem([]))
		for i in range(bits):
			for j in range(len(coeff) - 1, 0, -1):
				coeff[j] = coeff[j - 1] # shifting to the left -> we have little endian order!
			coeff[0] = GFElem([])
		coeff.reverse()
		return GFExtensionElem(coeff)

	def shiftRight(self, bits):
		coeff = []
		for i in range(len(self.coeff)):
			coeff.append(self.coeff[i])
		for i in range(bits):
			for j in range(len(coeff)):
				if (j < len(coeff) - 1):
					coeff[j] = coeff[j + 1] # shifting to the left -> we have little endian order!
					coeff[j + 1] = GFElem([])
		coeff.reverse()
		return GFExtensionElem(coeff)

	def isUnit(self):
		if (len(self.coeff) == 0):
			return False
		elif (self.coeff[0].isUnit()): #and len(self.coeff) == 1):
			for i in range(1, len(self.coeff)):
				if (not self.coeff[i].isZero()):
					return False
			return True
		else:
			return False

	def isZero(self):
		# for i in range(len(self.coeff)):
		# 	if (not self.coeff[i].isZero()):
		# 		return False
		for c in self.coeff:
			if (not c.isZero()):
				return False
		return True

	def getCoeff(self):
		return self.coeff[:]

	def bin(self, n):
		result = ""
		for c in self.coeff:
			subString = c.bin(fill = True, n = n)
			if (len(subString) < n):
				subString = ("0" * (n - len(subString))) + subString # CAW FLIPPED THIS
			result = result + subString[::-1]
		return result[::-1]

	def copy(self):
		c = []
		for i in range(len(self.coeff)):
			c.append(self.coeff[i].copy())
		c.reverse()
		return GFExtensionElem(c)

	def wt(self):
		count = 0
		for c in self.coeff:
			if (not c.isZero()):
				count = count + c.wt()
		return count

	def nterms(self):
		count = 0
		for c in self.coeff:
			if (not c.isZero()):
				count = count + 1
		return count

	def __eq__(self, other):
		if isinstance(other, GFExtensionElem):
			match = True
			length = max(len(self.coeff), len(other))
			for i in range(length):
				if not (self[i] == other[i]):
					match = False
					break
			return match

	def __len__(self):
		return len(self.coeff)

	def __getitem__(self, key):
		if (len(self.coeff) <= key):
			for i in range(key - len(self.coeff) + 1):
				self.coeff.append(GFElem([])) # 0
		return self.coeff[key]

	def __setitem__(self, key, value):
		if (len(self.coeff) <= key):
			for i in range(key - len(self.coeff) + 1):
				self.coeff.append(GFElem([])) # 0
		self.coeff[key] = value	

	def toString(self, trim = True):
		if (trim == False):
			return str(self)
		else:
			result = ""
			poly = ""
			self.coeff.reverse()
			exp = len(self.coeff) - 1
			zero = GFElem([])
			for x in range(len(self.coeff)):
				if not (self.coeff[x] == zero):
					if (str(exp) == "0"):
						poly = poly + self.coeff[x].toString(trim = trim) + " "
					elif (str(exp) == "1"):
						poly = poly + "(" + self.coeff[x].toString(trim = trim) + ")*y + "
					else:
						poly = poly + "(" + self.coeff[x].toString(trim = trim) + ")*y^" + str(exp) + " + "
				exp = exp - 1
			poly = poly[0:len(poly)]
			if (poly == ""):
				poly = "0"
			if (poly.endswith("+ ")):
				poly = poly[0:len(poly) - 2]
			self.coeff.reverse()
			return poly

	def __str__(self):
		poly = "["
		self.coeff.reverse()
		exp = len(self.coeff) - 1
		zero = GFElem([])
		for x in range(len(self.coeff)):
			if not (self.coeff[x] == zero):
				poly = poly + str(self.coeff[x]) + "y^" + str(exp) + " + "
			exp = exp - 1
		poly = poly[0:len(poly) - 3] + "]"
		if (poly == "]"):
			poly = "0"
		self.coeff.reverse()
		return poly

class GF:
	def __init__(self, base, exp, ip):
		self.base = base
		self.exp = exp
		self.ip = ip # must be instance of GFElem

		self.elems = []
		for i in range(2**exp):
			coeff = createBaseElem(2, i, exp)
			elem = GFElem(coeff)
			self.elems.append(elem)

	def g_add(self, x, y):
		sum = []
		index = 0
		for i in range(len(x) + len(y)):
			z = (x[i] + y[i]) % self.base
			sum.insert(0, z)
			index = index + 1
		result = GFElem(sum)
		return result

	def g_sub(self, x, y):
		diff = []
		index = 0
		for i in range(len(x) + len(y)):
			z = (x[i] - y[i]) % self.base
			diff.insert(0, z)
			index = index + 1
		return GFElem(diff)

	def g_mult(self, x, y, reduce = True):
		prod = []
		for i in range(len(x) + len(y)):
			prod.append(0)
		for i in range(len(x)):
			for j in range(len(y)): 
				pl = i + j
				tempProd = 0
				if (i < len(x) and j < len(y)):
					tempProd = x[i] * y[j]
				pl = i + j
				prod[pl] = prod[pl] + tempProd
		prod.reverse()

		# Reduce the elements not of the highest degree
		foundNonZero = False
		if (reduce == True):
			for i in range(len(prod)):
				prod[i] = prod[i] % self.base

		result = GFElem(prod)
		if (result.degree() >= self.exp and reduce == True):
			result = self.g_div(result, self.ip)[1] # div returns (Q,R), we want R(emainder)
		return result

	def reduce(self, x):
		result = x.copy()
		if (result.degree() >= self.exp):
			result = self.g_div(result, self.ip)[1] # div returns (Q,R), we want R(emainder)
		return result

	def g_div(self, N, D):
		if (D.degree() < 0):
			raise TypeError()
		if (N.degree() >= D.degree()): # numerator degree larger than divisor
			Q = GFElem([]) # q <- 0
			while (N.degree() >= D.degree()):
				divisor = D.shiftLeft(N.degree() - D.degree()) # find the divisor for the degrees - this is correct
				Q[N.degree() - D.degree()] = N[N.degree()] / divisor[divisor.degree()]
				divisor = divisor.scalarMult(Q[N.degree() - D.degree()]) 
				N = self.g_sub(N, divisor)
			R = N # left over numerator is the remainder...
			return (Q, R)
		else: # divisor degree larger than numerator
			Q = GFElem([]) # 0
			R = N
			return (Q, R)

	def inverse(self, x):
		g, x, y = self.EEA(x, self.ip)
		return self.reduce(x)

	def EEA(self, a, b):
		lastRem = a.copy()
		rem = b.copy()
		x = GFElem([0])
		lastx = GFElem([1])
		y = GFElem([1])
		lasty = GFElem([0])

		while not (rem.isZero()):
			# print(rem)
			tmp = rem.copy()
			quotient, rem = self.g_div(lastRem, rem)
			lastRem = tmp
			tmp = x.copy()
			x = self.g_sub(lastx, self.g_mult(quotient, x))
			lastx = tmp
			tmp = y.copy()
			y = self.g_sub(lasty, self.g_mult(quotient, y))
			lasty = tmp

		return lastRem, lastx, lasty

	def power(self, x, k):
		val = GFElem([1])
		exp = bin(k)[2:][::-1] # put the exponent in reverse order
		if k < 0:
			raise TypeError() 
		else:
			# for i in range(k):
			# 	val = self.g_mult(val, x)
			for i in range(len(exp)):
				if (exp[i] == '1'):
					val = self.g_mult(val, x)
				x = self.g_mult(x, x)
			return val

	def findOrder(self, g): # works for fields with base = 2 (characteristic)
		# p = []
		# gg = g
		# for j in range (self.exp):
		# 	p.append(gg % self.base)
		# 	gg = gg / self.base
		# p.reverse()

		coeff = createBaseElem(self.base, g, self.exp)
		gen = GFElem(coeff)
		curr = self.g_mult(gen, gen)
		cycled = False
		counter = 1
		coeffs = []
		coeffs.append(gen.getCoeff())
		while (not cycled):
			if (curr == gen or curr.isZero() or curr.getCoeff() in coeffs): 
				cycled = True
			else:
				counter = counter + 1
				coeffs.append(curr.getCoeff())
				curr = self.g_mult(gen, curr)
		return (gen, counter)

	def isGenerator(self, gen):
		# The field polynomial MUST be irreducible for this method to work, yes?
		order = self.getOrder() - 1
		factors = primeFactors(order)
		combs = []
		for i in factors:
			for j in factors:
				prod = i * j
				if (prod < order):
					combs.append(prod)
		if (gen.isUnit()):
			return False
		for i in factors:
			x = self.power(gen, i)
			if (x.isUnit() or x.isZero()):
				# print("Cyclic from " + str(i) + " factor: " + str(x))
				return False
		for i in combs:
			x = self.power(gen, i)
			if (x.isUnit() or x.isZero()):
				# print("Cyclic from " + str(i) + " combination: " + str(x))
				return False
		return True

	def findGeneratorsFast(self, earlyTerm = False, lb = -1, ub = -1):
		fieldOrder = self.getOrder()
		numGens = self.getBase() ** (self.getExp() - 1)
		gens = []
		count = 0
		mylb = 0
		myub = fieldOrder
		found = False
		if (lb != -1 and ub != -1):
			mylb = lb
			myub = ub
		for i in range(mylb, myub):
			coeff = createBaseElem(self.base, i, self.exp)
			gen = GFElem(coeff)
			print >> sys.stderr, "Trying: " + str(gen)
			if (self.isGenerator(gen)):
				gens.append(gen)
				if (earlyTerm):
					return gens
			count = count + 1
			# if (count > numGens and found == False): # yeah...
			# 	print >> sys.stderr, "Didn't encounter a generator within the logical threshold"
			# 	return gens
		return gens

	def findGenerators(self, earlyTerm = False, threshold = -1, lb = -1, ub = -1): # has to be within the first 100 elements...
		fieldOrder = self.getOrder()
		numGens = self.getBase() ** (self.getExp() - 1)
		gens = []
		count = 0
		mylb = 0
		myub = fieldOrder
		found = False
		if (lb != -1 and ub != -1):
			mylb = lb
			myub = ub
		print >> sys.stderr, "Finding generators " + str(lb) + " " + str(ub) + " " + str(mylb) + " " + str(myub)
		print >> sys.stderr, "Threshold: " + str(numGens)
		for i in range(mylb, myub):
			elemOrder = self.findOrder(i)
			print >> sys.stderr, str(i)
			if (elemOrder[1] == (fieldOrder - 1)): # it's a generator...
				gens.append(elemOrder[0]) 
				found = True
				if (earlyTerm): # early termination for polynomial generation code
					return gens
			count = count + 1
			if (threshold != -1 and count > threshold): # don't keep looking if we haven't found one yet...
				print >> sys.stderr, "Passed the specified threshold for finding generators in a base field"
				return gens
			if (count > numGens and found == False): # yeah...
				print >> sys.stderr, "Didn't encounter a generator within the logical threshold"
				return gens
		return gens

	def isRoot(self, x):
		accum = GFElem([]) # sum = 0 to start
		for i in range(len(self.ip)):
			print(i)
			if (self.ip.lsb(i) != 0):
				elem = self.power(x, i)
				accum = self.g_add(accum, elem)
				print(accum)
		return accum.isZero()

	def isExtensionRoot(self, x, eField): # we need the extension field to do arithmetic...
		accum = GFExtensionElem([]) # sum = 0 to start
		for i in range(len(self.ip)):
			print(i)
			if (self.ip.lsb(i) != 0):
				print(x)
				elem = eField.power(x, i)
				accum = eField.g_add(accum, elem)
				print(accum)
		return accum.isZero()

	def getElems(self):
		return self.elems

	def getIp(self):
		return self.ip

	def getOrder(self):
		return self.base ** self.exp

	def getBase(self):
		return self.base

	def getExp(self):
		return self.exp

	def hashcode(str):
		return self.ip.bin()

	def __str__(self):
		''' Format the field nicely for display.
		'''
		result = "GF(" + str(self.base) + "^" + str(self.exp) + ")"
		result = result + ", P(x) = " + str(self.ip)
		return result

class GFExtension:
	def __init__(self, baseField, extension, ip):
		self.baseField = baseField
		self.extension = extension
		self.ip = ip

	def g_add(self, x, y):
		sum = []
		for i in range(len(x) + len(y)):
			z = self.baseField.g_add(x[i], y[i])
			sum.insert(0, z)
		return GFExtensionElem(sum)

	def g_sub(self, x, y):
		diff = []
		for i in range(len(x) + len(y)):
			z = self.baseField.g_sub(x[i], y[i])
			diff.insert(0, z)
		return GFExtensionElem(diff)

	def g_mult(self, x, y):
		prod = []
		for i in range(len(x) + len(y)):
			prod.append(GFElem([]))
		for i in range(len(x)):
			for j in range(len(y)): 
				pl = i + j
				tempProd = GFElem([])
				if (i < len(x) and j < len(y)):
					tempProd = self.baseField.g_mult(x[i], y[j])
				pl = i + j
				prod[pl] = self.baseField.g_add(prod[pl], tempProd)
		prod.reverse()

		# Reduce the elements not of the highest degree
		foundNonZero = False
		for i in range(len(prod)):
			prod[i] = self.baseField.g_div(prod[i], self.baseField.getIp())[1]

		result = GFExtensionElem(prod)
		if (result.degree() >= self.extension):
			result = self.g_div(result, self.ip)[1] # div returns (Q,R), we want R(emainder)
		return result

	def g_div(self, N, D):
		if (D.degree() < 0): 
			raise TypeError() # cannot divide by zero.
		if (N.degree() >= D.degree()): # numerator degree larger than divisor
			Q = GFExtensionElem([GFElem([])]) # q <- 0
			Ns = []
			while (N.degree() >= D.degree()):
				divisor = D.shiftLeft(N.degree() - D.degree()) # find the divisor for the degrees - this is correct
				Q[N.degree() - D.degree()] = self.baseField.g_div(N[N.degree()], divisor[divisor.degree()])[0] # we want the quotient, not remainder
				divisor = divisor.scalarMult(Q[N.degree() - D.degree()], self.baseField)
				N = self.g_sub(N, divisor)
			R = N # left over numerator is the remainder...
			return (Q, R)
		else: # divisor degree larger than numerator
			Q = GFExtensionElem([]) # 0
			R = N
			return (Q, R)

	def findOrder(self, g): 
		p = []
		gg = g
		subOrder = self.baseField.getOrder()
		for j in range (self.extension):
			p.append(GFElem([], id = gg % subOrder, base = self.baseField.getBase(), exp = self.baseField.getExp()))
			gg = gg / subOrder
		p.reverse()
		gen = GFExtensionElem(p)
		curr = self.g_mult(gen, gen)
		cycled = False
		counter = 1
		coeffs = []
		n = self.baseField.getExp()
		coeffs.append(gen.bin(n))
		while (not cycled):
			bin = curr.bin(n)
			if (curr == gen or curr.isZero() or bin in coeffs): 
				cycled = True
			else:
				counter = counter + 1
				coeffs.append(bin)
				curr = self.g_mult(gen, curr)
		return (gen, counter)

	def isGenerator(self, gen):
		# The field polynomial MUST be irreducible for this method to work.
		order = self.getOrder() - 1
		factors = primeFactors(order)
		combs = []
		for i in factors:
			for j in factors:
				prod = i * j
				if (prod < order):
					combs.append(prod)
		if (gen.isUnit()):
			return False
		for i in factors:
			x = self.power(gen, i)
			if (x.isUnit() or x.isZero()):
				return False
		for i in combs:
			x = self.power(gen, i)
			if (x.isUnit() or x.isZero()):
				return False
		return True

	def findGeneratorsFast(self, earlyTerm = False, threshold = -1, lb = -1, ub = -1):
		fieldOrder = self.baseField.getOrder() ** self.extension
		numGens = self.baseField.getBase() ** ((self.baseField.getExp() * self.extension) - 1)
		gens = []

		# Counter to see if we passed the threshold... 
		# probably not useful to use this anymore
		count = 0
		found = False

		# Bounds for generator check
		mylb = 2 ** self.getN()
		myub = fieldOrder
		if (lb != -1 and ub != -1):
			mylb = lb
			myub = ub
		for i in range(mylb, myub): # start with at least one y coefficient (that's a requirement for a primitive element)
			try:
				p = []
				gg = i
				subOrder = self.baseField.getOrder()
				for j in range (self.extension):
					p.append(GFElem([], id = gg % subOrder, base = self.baseField.getBase(), exp = self.baseField.getExp()))
					gg = gg / subOrder
				p.reverse()
				gen = GFExtensionElem(p)
				print >> sys.stderr, "Trying extension fast generator: " + str(gen)
				if (self.isGenerator(gen)):
					print >> sys.stderr, "Generator."
					gens.append(gen)
					found = True
					if (earlyTerm):
						return gens
				else:
					print >> sys.stderr, "Not a generator."
				count = count + 1
				if (count > numGens and found == False): # yeah...
					print >> sys.stderr, "Didn't encounter a generator within the logical threshold"
					return gens
				if (threshold != -1 and count > threshold): # don'e keep looking if we haven't found one yet...
					print >> sys.stderr, "Passed the threshold..."
					return gens
			except Exception as e:
				print >> sys.stderr, "Error in GFExtension.findGenerators: " + str(e)
				# traceback.print_exc(limit=2,file=sys.stderr)
				pass
		return gens

	def findGenerators(self, earlyTerm = False, threshold = 2, lb = -1, ub = -1):
		fieldOrder = self.baseField.getOrder() ** self.extension
		numGens = self.baseField.getBase() ** ((self.baseField.getExp() * self.extension) - 1)
		gens = []

		# Counter to see if we passed the threshold... 
		# probably not useful to use this anymore
		count = 0

		# Bounds for generator check
		mylb = 2 ** self.getN()
		myub = fieldOrder
		if (lb != -1 and ub != -1):
			mylb = mylb + lb
			myub = mylb + ub
			if (myub > fieldOrder):
				myub = fieldOrder
		print >> sys.stderr, str(lb) + " " + str(ub) + " " + str(mylb) + " " + str(myub)
		for i in range(mylb, myub): # start with at least one y coefficient (that's a requirement for a primitive element)
			try:
				elemOrder = self.findOrder(i)
				print >> sys.stderr, "Trying extension generator index: " + str(i)
				if (elemOrder[1] == (fieldOrder - 1)): 
					print >> sys.stderr, "Generator."
					gens.append(elemOrder[0])
					if (earlyTerm):
						return gens
				else:
					print >> sys.stderr, "Not a generator."
				count = count + 1
				if (threshold != -1 and count > threshold): # don'e keep looking if we haven't found one yet...
					print >> sys.stderr, "Didn't encounter a generator within the logical threshold"
					return gens
			except Exception as e:
				# print >> sys.stderr, "Error in GFExtension.findGenerators: " + str(e)
				# traceback.print_exc(limit=2,file=sys.stderr)
				pass
		return gens

	def power(self, x, k):
		val = GFExtensionElem([GFElem([1])])
		exp = bin(k)[2:][::-1] # put the exponent in reverse order
		if k < 0:
			raise TypeError() 
		elif k == 0:
			return GFExtensionElem([GFElem([1])])
		else:
			# for i in range(k):
			# 	val = self.g_mult(e, x)
			for i in range(len(exp)):
				if (exp[i] == '1'):
					val = self.g_mult(val, x)
				x = self.g_mult(x, x)
			return val

	def getIP(self):
		return self.ip

	def getN(self):
		return self.baseField.getExp()

	def getM(self):
		return self.extension

	def getBaseField(self):
		return self.baseField

	def getExtension(self):
		return self.extension

	def getOrder(self):
		return (self.baseField.getBase() ** (self.baseField.getExp() + self.extension))

	def hashcode(str):
		return self.baseField.hashcode() + "-" + self.ip.bin(self.baseField.getExp())

	def __str__(self):
		''' Format the field nicely for display.
		'''
		result = str(self.extension) + " degree extension of: " + str(self.baseField) + ", using " + str(self.ip)
		return result

def test1():
    ip = GFElem([1,0,0,0,1,1,0,1,1]) #x^8 + x^4 + x^3 + x + 1
    field = GF(2, 8, ip)

    # Sample polynomials
    x = GFElem([1,1])
    y1 = GFElem([1,1])
    y2 = GFElem([1,0])
    y3 = GFElem([])
    y4 = GFElem([1,0,0,0,0,0,0,0])
    y5 = GFElem([1,0,0,0,0,0,0])
    v3 = GFElem([1,1])
    v4 = GFElem([1,0])
    v7 = GFElem([1,1,1])
    x10 = GFElem([1,1,1,1,1,1,1,1])
    y10 = GFElem([0,0,0,0,0,0,1,1])
    x11 = GFElem([0,1,0,1,0,0,1,1])
    y11 = GFElem([1,1,0,0,1,0,1,0])
    
    print("STARTING TEST CASES FOR GF(2^8)")

    # Sum test cases
    sum = field.g_add(x, y1)
    print(str(x) + " + " + str(y1) + " = " + str(sum))
    sum = field.g_add(x, y2)
    print(str(x) + " + " + str(y2) + " = " + str(sum))

    # Diff test cases
    diff = field.g_sub(x, y1)
    print(str(x) + " - " + str(y1) + " = " + str(diff))
    diff = field.g_sub(x, y2)
    print(str(x) + " - " + str(y2) + " = " + str(diff))

    # Mult test cases
    mult = field.g_mult(x, y1)
    print(str(x) + " * " + str(y1) + " = " + str(mult))
    mult = field.g_mult(v3, v7)
    print(str(v3) + " * " + str(v7) + " = " + str(mult))
    mult = field.g_mult(y5, y5)
    print(str(y5) + " * " + str(y5) + " = " + str(mult))
    mult = field.g_mult(y4, y4)
    print(str(y4) + " * " + str(y4) + " = " + str(mult))
    mult = field.g_mult(x10, y10)
    print(str(x10) + " * " + str(y10) + " = " + str(mult))
    mult = field.g_mult(x11, y11)
    print(str(x11) + " * " + str(y11) + " = " + str(mult))

    # Division cases
    div = field.g_div(x, ip)
    print(str(x) + " / " + str(ip) + " = " + str(div[0]) + "," + str(div[1]))

    print("STARTING TEST CASES FOR GF(2^2)")

    ip = GFElem([1, 1, 1]) #x^2 + x + 1
    field = GF(2, 2, ip)
    print(field)
    x0 = GFElem([0, 0])
    x1 = GFElem([0, 1])
    x2 = GFElem([1, 0])
    x3 = GFElem([1, 1])
    mult = field.g_mult(x0, x1)
    print(str(x0) + " * " + str(x1) + " = " + str(mult))
    mult = field.g_mult(x0, x2)
    print(str(x0) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x0, x3)
    print(str(x0) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x1, x2)
    print(str(x1) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x1, x3)
    print(str(x1) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x2, x3)
    print(str(x2) + " * " + str(x3) + " = " + str(mult))

    print("STARTING TEST CASES FOR GF(4^2)")

    ip = GFElem([1, 1, 1]) #x^2 + x + 1
    field = GF(4, 2, ip)
    print(field)
    x0 = GFElem([0, 0])
    x1 = GFElem([0, 1])
    x2 = GFElem([1, 0])
    x3 = GFElem([1, 1])
    mult = field.g_mult(x0, x1)
    print(str(x0) + " * " + str(x1) + " = " + str(mult))
    mult = field.g_mult(x0, x2)
    print(str(x0) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x0, x3)
    print(str(x0) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x1, x2)
    print(str(x1) + " * " + str(x2) + " = " + str(mult))
    mult = field.g_mult(x1, x3)
    print(str(x1) + " * " + str(x3) + " = " + str(mult))
    mult = field.g_mult(x2, x3)
    print(str(x2) + " * " + str(x3) + " = " + str(mult))

def test2():
	# neIp = GFElem([1,0,0,0,1,1,0,1,1]) #x^8 + x^4 + x^3 + x + 1
	# neField = GF(2, 8, neIp) # GF(2^8)
	# ip = GFElem([1,0,0,1,1]) # x^4 + x + 1
	# smallField = GF(2, 4, ip)
	# eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])]) # x^2 + x + 1100
	# eField = GFExtension(smallField, 2, eIp) # GF((2^4)^2)
	# gens = eField.findGenerators()
	# for g in gens:
	# 	print(g)

	

	ip = GFElem([1,0,0,1,1]) # x^4 + x + 1
	smallField = GF(2, 4, ip)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])]) # x^2 + x + 1100
	eField = GFExtension(smallField, 2, eIp) # GF((2^4)^2)
	gens = eField.findGenerators()
	for g in gens:
		print(g)

	# x = GFElem([1,1,0,1,0,0,0,0])
	# xInv = neField.inverse(x)
	# print("x and inverse")
	# print(x)
	# print(xInv)
	# mult = neField.g_mult(x, xInv)
	# print(mult)
	# y = GFElem([1,0,0,0,1,0,0])
	# # mult = neField.g_mult(x, y)
	# # print(str(x) + " * " + str(y) + " = " + str(mult))

	# x = GFExtensionElem([GFElem([1,1,1,0]), GFElem([1,0,1,0])])
	# y = GFExtensionElem([GFElem([1,1,1]), GFElem([1,0,0,0])])
	# # mult = eField.g_mult(x, y)
	# # print(str(x) + " * " + str(y) + " = " + str(mult))

	# N = GFExtensionElem([GFElem([1,1,0,0]), GFElem([1,0,1,0]), GFElem([1,1,1,1])])
	# D = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])])
	# # y = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])])
	# # print(y.shiftLeft(1))
	# # print(y.shiftLeft(2))

	# ip = GFElem([1,1,1])
	# smallField = GF(2,2,ip)
	# eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,0])])
	# eField = GFExtension(smallField, 2, eIp)
	# print(eField)
	# gens = eField.findGenerators()
	# print("Displaying generators:")
	# for g in gens:
	# 	print(g)

def testIsRoot():
	ip = GFElem([1,0,0,1,1]) # 4 1 0
	bigField = GF(2,4,ip)
	neIp = GFElem([1,1,1]) # 2 1 0
	neField = GF(2,2,neIp)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1, 0])]) # 2 1 {1}0
	eField = GFExtension(neField, 2, eIp)
	print(bigField)
	print(eField)
	print(bigField.isRoot(GFElem([1,0])))
	print(bigField.isExtensionRoot(GFExtensionElem([GFElem([1]),GFElem([0])]), eField))

def findBigRoot():
	ip = GFElem([]) # 4 1 0
	bigField = GF(2,16,ip)
	neIp = GFElem([1,0,0,0,1,1,1,0,1]) # 8 4 3 2 0
	neField = GF(2,8,neIp)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1, 0])]) # 2 1 {1}0
	eField = GFExtension(neField, 2, eIp)

def testPower():
	''' Simple test for the efficiency of the new multiplication algorithm.
	'''
	ip = GFElem([1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1]) # x^4 + x + 1
	smallField = GF(2, 16, ip)
	x = GFElem([1,0,1])

	start = time.time()
	val = smallField.g_mult(x, x)
	for i in range(100):
		val = smallField.g_mult(x, val)
	end = time.time()

	start2 = time.time()
	val2 = smallField.power(x, 102)
	end2 = time.time()

	print(str(x) + " cubed = " + str(val) + ", " + str(end - start))
	print(str(x) + " cubed = " + str(val2) + ", " + str(end2 - start2))

	# Now test for extension fields...
	ip = GFElem([1,0,0,1,1]) # x^4 + x + 1
	smallField = GF(2, 4, ip)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,1,0,0])]) # x^2 + x + 1100
	eField = GFExtension(smallField, 2, eIp) # GF((2^4)^2)

	x = GFExtensionElem([GFElem([1,0,1]),GFElem([1,1,1])])

	start = time.time()
	val = eField.g_mult(x, x)
	for i in range(100):
		val = eField.g_mult(x, val)
	end = time.time()

	start2 = time.time()
	val2 = eField.power(x, 102)
	end2 = time.time()

	print(str(x) + " cubed = " + str(val) + ", " + str(end - start))
	print(str(x) + " cubed = " + str(val2) + ", " + str(end2 - start2))

def testIsGenerator():
	# ip = GFElem([1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1]) # x^4 + x + 1
	# smallField = GF(2, 16, ip)
	ip = GFElem([1,0,0,0,1,1,0,1,1]) #x^8 + x^4 + x^3 + x + 1
	field = GF(2, 8, ip) # GF(2^8)
	gen  = GFElem([1,1])
	print(field.isGenerator(gen))
	gen  = GFElem([1,0])
	print(field.isGenerator(gen))
	print(len(field.findGenerators()))

def test3():
	# test = GFExtensionElem([GFElem([1]), GFElem([1])])
	# print(test.isUnit())
	# print(test.isZero())
	# zero = GFExtensionElem([GFElem([])])
	# print(zero.isUnit())
	# print(zero.isZero())
	# unit = GFExtensionElem([GFElem([1])])
	# print(unit.isUnit())
	# print(unit.isZero())
	# return

	neIp = GFElem([1,1,1]) # 8 4 3 2 0
	neField = GF(2,2,neIp)
	eIp = GFExtensionElem([GFElem([1]), GFElem([1]), GFElem([1,0])]) 
	eField = GFExtension(neField, 2, eIp)
	x = GFExtensionElem([GFElem([1,0]), GFElem([1])]) # 1 0 0 1
	print(neField)
	print(eField)
	print(x)
	prod = eField.g_mult(x, x)
	print(prod)
	print("")

	ip = GFElem([1,0,0,1,1]) # 4 1 0
	field = GF(2, 4, ip)
	x = GFElem([1,1,0,1])
	print(field)
	print(x)
	prod = field.g_mult(x, x)
	print(prod)

# Main entry point - just runs one of the tests if this file
# is ever started directly
def main():
	test3()

if __name__ == "__main__":
	main()
