/*
 * File: galois.c
 * Author: Christopher A. Wood, caw4567@rit.edu (www.christopher-wood.com)
 * Description: Implementation of GF(2^8) and GF(2^16) arithmetic
 * 	in standard polynomial basis and towers of normal bases. The ability
 * 	to use polynomial bases for the tower field representation is also
 * 	supported, pending the inclusion of an appropriate pair of basis
 * 	transformation matrices T and TINV (see comments). Furthermore,
 * 	the "standard" field (e.g. the AES field) irreducible polynomials
 * 	are fixed here. They may be easily changed, however. 
 */

#include "galois.h"

// Toggle these macros to set what's being tested in this module
// Read the comments for specifics on each tower-field construction
// and the relevant arithmetic.
#define VERSION_GF256_1    1
#define VERSION_GF256_2    0 
#define VERSION_GF6K_1     0

// Standard definitions (for sbox16.c)
#define PX_16  0x002B   
#define FPX_16 0x1002B  
#define PX_8   0x77    
#define FPX_8  0x177    

/** 
 * Standard test for GF(2^16) inverse calculation. 
 * -> T(v)   = v^16 + v^5 + v^3 + v + 1
 * -> Sigma  = v
 * -> Pi     = vw^4 + v^2 + v
 * -> Lambda = (v^2 + v)wx^16 + vw^4x
 * -> Bases  = [v, v^2], [w, w^4], [x, x^16], [y, y^256]
 */
#if VERSION_GF6K_1

#define GF4_NORMAL_BASIS
#define GF16_NORMAL_BASIS
#define GF256_NORMAL_BASIS
#define GF6K_NORMAL_BASIS

#define PX_16  0x002B   
#define FPX_16 0x1002B  // include the MSB (1) V^16 + V^5 + V^3 + V + 1
#define PX_8   0x77     // NOT USED
#define FPX_8  0x177    // NOT USED

static uint16_t TINV[16] = 
{
	0xa01a, 0x599f, 0x9eb4, 0x52e9, 0xc3c4, 0xb5a4, 0x018a, 0x4a64,
	0xce25, 0x98ba, 0x8dfa, 0xa16b, 0x3fd1, 0x8649, 0x7871, 0xffff
};

static uint16_t T[16] = 
{
	0xd31e, 0xfcf5, 0x5f95, 0x63ea, 0x26b9, 0xc683, 0xc6db, 0x6621,
	0x957c, 0xc328, 0x0ffa, 0xfdd6, 0xc75f, 0xd145, 0x3602, 0xd735
};

// Define the coefficient values
#define SIGMA   0x1   
#define PI      0x7   
#define LAMBDA  0x34  
#endif

/** 
 * Test for GF(2^8) inverse calculation using AES polynomial.
 * 	This uses an alternate isomorphic mapping from GF(2^8) to GF(((2^2)^2)^2).
 *  See the thesis document for details.
 * -> S(v)   = v^8 + v^4 + v^3 + v + 1
 * -> Sigma  = v^2
 * -> Pi     = vw
 * -> Bases  = [v, v^2], [w, w^4], [x^16, x]
 */
#if VERSION_GF256_1

#define GF4_NORMAL_BASIS
#define GF16_NORMAL_BASIS
#define GF256_NORMAL_BASIS
#define GF6K_NORMAL_BASIS

#define PX_16  0x002B   // NOT USED
#define FPX_16 0x1002B  // NOT USED
#define PX_8   0x1B     
#define FPX_8  0x11B    // include the MSB (1)

static uint8_t TINV[8] = 
{
	0x98, 0xf3, 0xf2, 0x48, 0x09, 0x81, 0xa9, 0xff
};
static uint8_t T[8] = 
{
	0x64, 0x78, 0x6e, 0x8c, 0x68, 0x29, 0xde, 0x60
};

// Define the coefficient values
#define SIGMA   0x2   
#define PI      0x1
#define LAMBDA  0x34   // NOT USED   
#endif

/** 
 * Test for GF(2^8) inverse calculation using new AES polynomial 
 * -> S(v)   = v^8 + v^6 + v^5 + v^4 + v^2 + v + 1
 * -> Sigma  = v
 * -> Pi     = (v + 1)w^4 + vw
 * -> Bases  = [1, v], [w, w^4], [x^16, x]
 */
#if VERSION_GF256_2

#define GF4_POLYNOMIAL_BASIS
#define GF16_NORMAL_BASIS
#define GF256_NORMAL_BASIS
#define GF6K_NORMAL_BASIS

#define PX_16  0x002B   // NOT USED
#define FPX_16 0x1002B  // NOT USED
#define PX_8  0x77
#define FPX_8 0x177     // include the MSB (1)

static uint8_t TINV[8] =
{
	0x3a, 0x13, 0x64, 0x01, 0x2b, 0x4a, 0xea, 0x55
};
static uint8_t T[8] = 
{
	0xaf, 0xb5, 0xa9, 0x98, 0x79, 0x3c, 0xc8, 0x10
};

// Define the coefficient values
#define SIGMA   0x2   
#define PI      0xD
#define LAMBDA  0x0 // UNUSED
#endif

uint8_t g8_add(uint8_t x, uint8_t y) 
{
	return x ^ y;
}

uint8_t g8_mul(uint8_t x, uint8_t y) 
{
	int i;
	uint8_t sum = 0;
	int msbSet;
	for (i = 0; i < 8; ++i) 
	{
		if (1 & y) sum ^= x;
		msbSet = x & MSB_8;
		x <<= 1;
		if (msbSet) x ^= FPX_8 & 0xFF;
		y >>= 1;
	}
	return sum;
}

QR g8_div(uint16_t ai, uint8_t b) 
{
	uint8_t a = ai;
	int msb = MSB_8;
	int d = 0;
	QR result = {0, 0};
	while (b > 0 && !(b & MSB_8)) 
	{
		++d;
		b <<= 1;
	}
	if (ai & HMSB_8) 
	{
		result.q ^= 1 << (d+1);
		a ^= b << 1;
	}
	for (;d > -1; d--) 
	{
		if ((a & msb) && (b & msb)) 
		{
			result.q ^= 1 << d;
			a ^= b;
		}
		msb >>= 1;
		b >>=1;
	}
	result.r = a;
	return result;
}

uint8_t g8_inv(uint8_t x) 
{
	if (x == 0) return 0;
	if (x == 1) return 1;

	uint8_t r0 = PX_8; // rem[i - 2]
	uint8_t r1 = x;    // rem[i - 1]
	uint8_t a0 = 0;    // aux[i - 2]
	uint8_t a1 = 1;    // aux[i - 1]
	uint8_t tmp;
	QR qr;

	int firstRun = 0;
	while (r1 > 0)
	{
		if (firstRun != 0) qr = g8_div(r0, r1);
		else
		{
			qr = g8_div(FPX_8, r1); 
			firstRun++;
		}
		r0 = r1; r1 = qr.r;
		tmp = a0; a0 = a1;
		a1 = g8_add(tmp, g8_mul(qr.q, a1));
	}

	return a0;
}

uint16_t g16_add(uint16_t x, uint16_t y)
{
	return x ^ y;
}

uint16_t g16_sub(uint16_t x, uint16_t y)
{
	return x ^ y;
}

uint16_t g16_mul(uint16_t x, uint16_t y)
{
	uint16_t accum = 0;
	uint16_t msb = 0;
	uint16_t i;
	for (i = 0; i < 16; i++) 
	{
		if (y & LSB) accum ^= x;
		msb = (x & MSB_16); // fetch the MSB
		x <<= 1;
		if (msb) x ^= PX_16;
		y >>= 1;
	}
	return accum;
}

/**
 * Polynomial division in GF(2^16).
 */
QR g16_div(uint32_t ai, uint16_t b)
{
	uint16_t a = (uint16_t)ai;
	int msb = MSB_16;
	int d = 0;
	QR result = {0, 0};

	// Align the denominator with the numerator
	while (b > 0 && !(b & MSB_16)) {
		++d;
		b <<= 1;
	}

	// If the polynomial MSB is set (17th bit), increment 
	// the quotient and reduce the numerator.
	if (ai & HMSB_16) {
		result.q ^= 1 << (d+1);
		a ^= b << 1;
	}

	for (;d > -1; d--) {
		if ((a & msb) && (b & msb)) {
			result.q ^= 1 << d;
			a ^= b;
		}
		msb >>= 1;
		b >>=1;
	}

	result.r = a;
	return result;
}

/**
 * Modular inverse in GF(2^16) using the EEA algorithm.
 */
uint16_t g16_inv(uint16_t x)
{
	// Trivial special cases.
	if (x == 0) return 0;
	if (x == 1) return 1;

	uint16_t r0 = PX_16; // rem[i - 2]
	uint16_t r1 = x;     // rem[i - 1]
	uint16_t a0 = 0;     // aux[i - 2]
	uint16_t a1 = 1;     // aux[i - 1]
	uint16_t tmp;
	QR qr;

	int firstRun = 0;
	while (r1 > 0)
	{
		if (firstRun != 0) qr = g16_div(r0, r1);
		else 
		{
			qr = g16_div(FPX_16, r1); 
			firstRun++;
		}
		r0 = r1; r1 = qr.r;
		tmp = a0; a0 = a1;
		a1 = g16_add(tmp, g16_mul(qr.q, a1));
	}

	return a0;
}

uint16_t g16_change_basis(uint16_t x, uint16_t* M)
{
	int32_t i;
	uint16_t y = 0;

	for (i = 15; i >= 0; i--) 
	{
		if (x & 1) y ^= M[i];
		x >>= 1;
	}

	return y;
}

uint8_t g8_change_basis(uint8_t x, uint8_t* M)
{
	int32_t i;
	uint8_t y = 0;

	for (i = 7; i >= 0; i--) 
	{
		if (x & 1) y ^= M[i];
		x >>= 1;
	}

	return y;
}

uint8_t g22_mul(uint8_t x, uint8_t y)
{
#ifdef GF4_NORMAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0x2) >> 1; 
	b = (x & 0x1);
	c = (y & 0x2) >> 1; 
	d = (y & 0x1);
	e = (a ^ b) & (c ^ d);
	p = (a & c) ^ e;
	q = (b & d) ^ e;
	return ((p << 1) | q);
#else // POLYNOMIAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0x2) >> 1; 
	b = (x & 0x1);
	c = (y & 0x2) >> 1; 
	d = (y & 0x1);
	e = (a ^ b) & (c ^ d);
	p = (b & d) ^ e;
	q = (a & c) ^ (b & d);
	return ((p << 1) | q);
#endif
}

uint8_t g22_sq(uint8_t x) 
{
	return g22_mul(x, x);
}

uint8_t g222_mul(uint8_t x, uint8_t y) 
{
#ifdef GF16_NORMAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xC) >> 2; 
	b = (x & 0x3);
	c = (y & 0xC) >> 2; 
	d = (y & 0x3);
	e = g22_mul(a ^ b, c ^ d);
	e = g22_mul(e, SIGMA);
	p = g22_mul(a, c) ^ e;
	q = g22_mul(b, d) ^ e;
	return ((p <<2 ) | q);
#else // POLYNOMIAL_BASIS
	uint8_t a, b, c, d, e, f, p, q;
	a = (x & 0xC) >> 2; 
	b = (x & 0xC);
	c = (y & 0xC) >> 2; 
	d = (y & 0xC);
	e = gf22_mul(a ^ b, c ^ d);
	f = gf22_mul(a ^ c, SIGMA);
	p = gf22_mul(b, d) ^ e;
	q = gf22_mul(b, d) ^ f;
	return ((p << 2) | q);
#endif
}

uint8_t g222_sq(uint8_t x) 
{
	return g222_mul(x, x);
}

uint8_t g222_inv(uint8_t x) 
{
#ifdef GF16_NORMAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xC) >> 2; 
	b = (x & 0x3);
	c = g22_mul(g22_sq(a ^ b), SIGMA);
	d = g22_mul(a, b);
	e = g22_sq(c ^ d);
	p = g22_mul(e, b);
	q = g22_mul(e, a);
	return ((p << 2) | q);
#else // POLYNOMIAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xC) >> 2; 
	b = (x & 0x3); 
	c = g22_mul(g22_sq(a ^ b), SIGMA); 
	d = a ^ b; 
	e = g22_sq(c ^ g22_mul(d, b)); // inverse
	p = g22_mul(e, a);
	q = g22_mul(e, d);
	return ((p << 2) | q);
#endif
}

uint8_t g2222_mul(uint8_t x, uint8_t y) 
{
#ifdef GF256_NORMAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xF0) >> 4; 
	b = (x & 0x0F);
	c = (y & 0xF0) >> 4; 
	d = (y & 0x0F);
	e = g222_mul(a ^ b, c ^ d);
	e = g222_mul(e, PI); 
	p = g222_mul(a, c) ^ e;
	q = g222_mul(b, d) ^ e;
	return ((p << 4) | q);
#else // POLYNOMIAL_BASIS
	uint8_t a, b, c, d, e, f, p, q;
	a = (x & 0xF0) >> 2; 
	b = (x & 0x0F);
	c = (y & 0xF0) >> 2; 
	d = (y & 0x0F);
	e = gf222_mul(a ^ b, c ^ d);
	f = gf222_mul(a ^ c, PI);
	p = gf222_mul(b, d) ^ e;
	q = gf222_mul(b, d) ^ f;
	return ((p << 4) | q);
#endif
}

uint8_t g2222_inv(uint8_t x) 
{
#ifdef GF256_NORMAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xF0) >> 4; 
	b = (x & 0x0F);
	c = g222_mul(g222_mul(a ^ b, a ^ b), PI); 
	d = g222_mul(a, b);
	e = g222_inv(c ^ d);
	p = g222_mul(e, b);
	q = g222_mul(e, a);
	return ((p << 4) | q);
#else // POLYNOMIAL_BASIS
	uint8_t a, b, c, d, e, p, q;
	a = (x & 0xC) >> 2; 
	b = (x & 0x3); 
	c = g222_mul(g222_mul(a ^ b, a ^ b), SIGMA); 
	d = a ^ b; 
	e = g222_inv(c ^ g222_mul(d, b)); // inverse
	p = g222_mul(e, a);
	q = g222_mul(e, d);
	return ((p << 4) | q);
#endif
}

uint16_t g22222_inv(uint16_t x) 
{
#ifdef GF6K_NORMAL_BASIS
	uint16_t a, b, c, d, e, p, q;
	a = (x & 0xFF00) >> 8; 
	b = (x & 0x00FF);
	c = g2222_mul(g2222_mul(a ^ b, a ^ b), LAMBDA); 
	d = g2222_mul(a, b);
	e = g2222_inv(c ^ d);
	p = g2222_mul(e, b);
	q = g2222_mul(e, a);
	return ((p << 8) | q);
#else // POLYNOMIAL_BASIS
	uint16_t a, b, c, d, e, p, q;
	a = (x & 0xFF00) >> 8; 
	b = (x & 0x00FF);
	c = g2222_mul(g2222_mul(a ^ b, a ^ b), LAMBDA); 
	d = a ^ b; 
	e = g2222_inv(c ^ g2222_mul(d, b)); // inverse
	p = g2222_mul(e, a);
	q = g2222_mul(e, d);
	return ((p << 8) | q);
#endif
}

//// MAIN ENTRY POINT ////
// Toggling the macros below plugs in the appropriate base change matrices and 
// coefficients to make the arithmetic work out in the respective basis.
// See the comments for each macro for a description of the tower field construction.

int main()
{
#if VERSION_GF256_1 | VERSION_GF256_2 // GF(2^8) test options
	uint8_t inv = 0;
	uint8_t cinv = 0;
	int x;
	for (x = 0; x <= 0xFF; x++)
	{
		x = (uint8_t)x;
		inv = g8_inv(x);
		cinv = g8_change_basis(x, TINV);
		cinv = g2222_inv(cinv);
		cinv = g8_change_basis(cinv, T);
		if (inv != cinv || (x != 0 && g8_mul(x, cinv) != 1))
		{
			printf("Failure.\n");
			return -1;
		}
	}
#if   VERSION_GF256_1
	printf("GF(((2^2)^2)^2) with basis [v,v^2],[w,w^4],[x^16,x] inverse success.\n");
#elif VERSION_GF256_2
	printf("GF(((2^2)^2)^2) with basis [1,v],[w,w^4],[x^16,x] inverse success.\n");
#endif 

#elif VERSION_GF6K_1 // GF(2^16) tests options
	uint16_t inv = 0;
	uint16_t cinv = 0;
	int x;
	for (x = 0; x <= 0xFFFF; x++)
	{
		x = (uint16_t)x;
		inv = g16_inv(x);
		cinv = g16_change_basis(x, TINV);
		cinv = g22222_inv(cinv);
		cinv = g16_change_basis(cinv, T);
		if (inv != cinv || (x != 0 && g16_mul(x, cinv) != 1))
		{
			printf("Failure.\n");
			return -1;
		}
	}
	printf("GF((((2^2)^2)^2)^2) with basis [v,v^2],[w,w^4],[x,x^16],[y,y^256] inverse success.\n");
#endif

	return 0;
}
