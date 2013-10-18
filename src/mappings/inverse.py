# File: inverse.py
# Author: Christopher A. Wood

import sys
sys.path.append("../CompositeFields/cfa/")
from galois import *
from galois_util import *

def usage():
	print >> sys.stderr, "Usage: python inverse.py <polynomial> <base> <extension>"

def main():
	if (len(sys.argv) != 4):
		usage()
		sys.exit(-1)
	else:
		# Read the command line arguments.
		poly = ""
		fbase = 0
		fextension = 0
		try:
			poly = sys.argv[1]
			fbase = int(sys.argv[2])
			fextension = int(sys.argv[3])
		except:
			print >> sys.stderr, "Invalid command line argument."
			usage()
			sys.exit(-1)

		# Create the field.
		ip = GFElem([], binString = poly)
		field = GF(fbase, fextension, ip)

		# Create the mapping.
		S = []
		order = fbase ** fextension
		for i in range(order):
			elem = GFElem([], id = i, base = fbase, exp = fextension)
			inverse = field.inverse(elem)

			# VERIFY
			if not (elem.isZero()):
				prod = field.g_mult(elem, inverse)
				assert(prod.isUnit())
			else:
				assert(inverse.isZero())

			S.append(inverse.bin())

		# Print the mapping
		out = ""
		for m in S:
			if (len(m) == 0):
				out = out + "0x0,"
			else:
				a = int(m,2)
				out = out + str(hex(a)) + ","
		print(out[0:len(out) - 1])

if __name__ == "__main__":
	main()
