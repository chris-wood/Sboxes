# File: galois_test.py
# Author: Christopher Wood

from galois import *

def main():
	y = GFElem([1,0,1,0,0,0,0,1,0,1,0,1,0])
	x = GFElem([1,0,1,0,0,1,0,0,0,0,1,1,1,1,1])
	z = GFElem([1,1,0,0,1,1,1,1,0,0,0])
	ip = GFElem([1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1])
	field = GF(2,16,ip)
	print(field)

	result = field.g_add(x,y)
	print(str(x) + " + " + str(y) + " = " + str(result))
	result = field.g_sub(x,y)
	print(str(x) + " - " + str(y) + " = " + str(result))
	result = field.g_mult(x,y)
	print(str(x) + " * " + str(y) + " = " + str(result))
	(q,r) = field.g_div(x,y)
	print(str(x) + " / " + str(y) + " = " + str(q) + ", " + str(r))
	print("----")
	result = field.inverse(x)
	print(str(x) + " ^{-1} = " + str(result))
	mult = field.g_mult(result, x)
	print("is it the inverse? " + str(mult))

if __name__ == "__main__":
	main()