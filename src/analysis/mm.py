# File: mm.py
# Author: Christopher Wood
# Description: Construction code for MM Boolean function
# NOTE: all polynomials are stored as little endian lists of coefficients

import sys
import pickle
sys.path.append("../CompositeFields/cfa/")
from galois import *
from galois_util import *

class LCode:
	def __init__(self, basis):
		self.u = len(basis[0])
		self.m = len(basis)
		self.basis = basis[:] # make a copy, don't store a reference

def xor(elems):
	result = []
	for i in range(len(elems[0])):
		val = 0
		for j in range(len(elems)):
			val = val ^ elems[j][i]
		result.append(val)
	return result

def apn54():
	# x -> x^3 in GF(2^5)
	ip = GFElem([1,0,0,0,1,1]) # x^5 + x + 1
	field = GF(2, 5, ip)
	powers = []
	for i in range(2**5):
		coeff = createBaseElem(2, i, 5)
		elem = GFElem(coeff)
		# e2 = GFElem(coeff)
		elem = field.power(elem, 3)
		# elem = field.g_mult(elem, e2)
		# elem = field.g_mult(elem, e2)
		bin = elem.bin(fill = True, n = 5)
		powers.append(bin) # [0:len(bin) - 1]
	return powers

def dMatrix(code, field, gen):

	print(field)
	print(gen)

	rows = (2 ** code.m - 1)

	powers = []
	for i in range(rows):
		alpha = field.power(gen, i)
		powers.append(alpha)

	D = {}
	for i in range(rows):
		entry = GFElem([])
		for j in range(4):
			elems = []
			for k in range(4):
				if (powers[(i + j) % rows][k] == 1):
					elems.append(code.basis[k])
			D[(i,j)] = xor(elems)

	# Debug...
	for i in range(rows):
		for j in range(4):
			print(str(D[(i,j)]) + "   "),
		print("")

	# Hard-code this as an example...
	L = [(2,35),(13,34)] # this is constructed based on the code - see the full paper
	print(apn54())

def main():

	mapping = apn54()
	for m in mapping:
		print(hex(int(m, 2)) + ","),

	code = LCode([GFElem([1,1,0,1,0,0,0]), GFElem([1,0,1,0,1,0,0]), GFElem([0,1,1,0,0,1,0]), GFElem([1,1,1,0,0,0,1])])
	ip = GFElem([1,0,0,1,1])
	field = GF(2,4,ip)
	gens = field.findGenerators(earlyTerm = True)
	dMatrix(code, field, gens[0])

if __name__ == "__main__":
	main()

