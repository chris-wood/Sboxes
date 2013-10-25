import sys

f = open(sys.argv[1], 'r')
count = -1
elems = []
for l in f:
	l = l.strip()
	if count >= 0:
		count = count + 1
	if l.startswith("*** Constant Inverse"):
		count = 0
	if l.startswith("*** Performing interpolation"):
		count = -1
	if count >= 2:
		elems.append(l)

for i in range(0, len(elems), 2):
	print >> sys.stderr, "e    = " + str(elems[i])
	print(elems[i])
	print >> sys.stderr, "e^-1 = " + str(elems[i + 1])
	print(elems[i+1])
