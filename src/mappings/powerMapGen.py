#####  #!/usr/bin/env sage -python

import sys
from sage.all import *

goon = True

#n = 8
n = int(sys.argv[1])

#x = GF(2,'a').polynomial()
#k = GF(2**n, name='x', modulus=x**8 + x**4 + x**3 + x + 1)
k = GF(2**n, 'x')

# from fractions import *
bijections = []
order = (2**n) - 1
for i in range(1, 2**n):
	if (gcd(i, order) == 1):
		bijections.append(i)
print(len(bijections))
# bijections = [2**n - 2]

if goon:
	for d in bijections:
		f = open('powermap_' + str(d) + '.out', 'w')
		f.write(str(k.polynomial()) + "\n")
		f.write(str(n) + " " + str(d) + "\n")
		print >> sys.stderr, "Generating x^" + str(d)
		bucket = []
		for e in k:
			if (e == 0):
				bucket.append((0, '0x0,'))
				# f.write("0\n")
				# f.write("0\n")	
			else:
				# Build input element hex string
				result = str(e)
				data = result.split(" + ")
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
				index = int(bsi, 2)

				# Build output element hex string
				val = e**d
				# f.write(str(e) + "\n")
				# f.write(str(val) + "\n")
				result = str(val)
				data = result.split(" + ")
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
				if index == 2:
					print(str(hex(int(bsi, 2))))
				bucket.append((index, str(hex(int(bsi, 2))) + ","))

		swapped = True
		while swapped:
			swapped = False
			for i in range(len(bucket)):
				for j in range(i + 1, len(bucket)):
					if (i != j):
						if (bucket[i][0] > bucket[j][0]):
							tmp = bucket[i]
							bucket[i] = bucket[j]
							bucket[j] = tmp
							swapped = True

		orderedBucket = []
		for i in range(len(bucket)):
			orderedBucket.append(bucket[i][1])

		mapping = ''.join(orderedBucket)
		f.write(mapping)
