# File: powerMapParse.py

import sys

def parsePoly(n, p):
	data = p.split(" + ")
	print(p)
	bs = []
	for i in range(n):
		bs.append(0)
	for b in data:
		if (b == '1'):
			bs[n - 1] = 1
		elif (b == '0'):
			bs[n - 1] = 0
		elif (b == 'x'):
			bs[n - 2] = 1
		else:
			bs[n - int(b[2:]) - 1] = 1
	bsi = ""
	for b in bs:
		bsi = bsi + str(b)
	val = int(bsi, 2)
	return val, str(hex(int(bsi, 2)))

def parsePowerMap(fname):
	f = open(fname, 'r')
	poly = f.readline().strip()
	header = f.readline().strip().split(" ")
	n, d = int(header[0]), int(header[1])

	# Read in the map, do the conversion, etc etc
	domain = []
	e = f.readline().strip()
	bucket = []
	while (len(e) > 0):
		ie = e
		e = f.readline().strip()
		oe = e
		iei, ieh = parsePoly(n, ie)
		oei, oeh = parsePoly(n, oe)
		bucket.append((iei, oeh + ","))

		# Advance...
		e = f.readline().strip()

	# Sort the domain map by the input values (their integer indices)
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
	of = open(fname + ".hex", 'w')
	of.write(mapping)

parsePowerMap(sys.argv[1])