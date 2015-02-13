/**
 * File: sbox16_unopt.c
 * Author: Christopher A. Wood, caw4567@rit.edu
 * Description: Software implementation of our unoptimized 16-bit S-box.
 * 	See the thesis for details on its construction.
 */

#include "galois.h"

static uint16_t A[16] =
{
	0x797b, 0x7c85, 0x9378, 0x151, 0x2312, 0x82f, 0x3f35, 0xe57e, 
	0x29d, 0x7e12, 0xdc62, 0xadbb, 0xced3, 0x87a0, 0xe900, 0x2d9c
};

static uint16_t AINV[16] =
{
	0x6a5e, 0xc863, 0x3b62, 0xec10, 0x3931, 0xb56e, 0xd1e7, 0xa06c,
	0x585f, 0x230c, 0xf6e0, 0x5557, 0x577e, 0x4d26, 0x17be, 0xc637
};

// Constant from our new S-box
#define C 0x45B7

uint16_t forward(uint16_t x)
{
	uint16_t inv = g16_inv(x);
	inv = g16_change_basis(inv, A);
	return inv ^ C;
}

uint16_t inverse(uint16_t x)
{
	uint16_t inv = x ^ C;
	inv = g16_change_basis(inv, AINV);
	return g16_inv(inv);
}

int main()
{
	uint16_t x, inv;
	uint16_t S[0x10000];
	uint16_t SINV[0x10000];
	int i, j;

	// Compute the S-box values
	for (i = 0; i <= 0xFFFF; i++)
	{
		x = (uint16_t)i;
		S[i] = forward(x);
		SINV[i] = inverse(x);
		if (inverse(forward(x)) != x)
		{
			printf("FAILURE with 0x%x\n", x);
			return -1;
		}
	}

	// Display the tables for each.
	printf("uint16_t S[65536] = {\n");
	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < 256; j++)
		{
			if (j == 255) printf("%x\n", forward((256 * i) + j));
			else printf("%x, ", forward((256 * i) + j));
		}
	}
	printf("};\n");
	printf("uint16_t SI[65536] = {\n");
	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < 256; j++)
		{
			if (j == 255) printf("%x\n", inverse((256 * i) + j));
			else printf("%x, ", inverse((256 * i) + j));
		}
	}
	printf("};\n");
}
