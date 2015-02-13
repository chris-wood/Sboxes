#include "galois.h"

static uint8_t A[8] =
{
	0x4d,
	0x58,
	0x20,
	0x4,
	0xd7,
	0xb2,
	0x71,
	0xf0
};

static uint8_t AINV[8] =
{
	0xcd,
	0xf2,
	0x20,
	0x1e,
	0xac,
	0x10,
	0xf7,
	0xce
};

// Include in galois.c
// #define PX_8  0x77
// #define FPX_8 0x177

// constant from our new S-box
#define C 0x08

uint8_t forward(uint8_t x)
{
	uint8_t inv = g8_inv(x);
	inv = g8_change_basis(inv, A);
	return inv ^ C;
}

uint8_t inverse(uint8_t x)
{
	uint8_t inv = x ^ C;
	inv = g8_change_basis(inv, AINV);
	return g8_inv(inv);
}

int main()
{
	uint8_t x, inv;

	uint8_t S[0x100];
	uint8_t SINV[0x100];

	int i, j;
	for (i = 0; i <= 0xFF; i++)
	{
		x = (uint8_t)i;
		S[i] = forward(x);
		SINV[i] = inverse(x);
		// printf("%x %x %x\n", x, forward(x), inverse(forward(x)));
		printf("0x%x -> 0x%x\n", x, forward(x));
		if (inverse(forward(x)) != x)
		{
			printf("FAILURE with 0x%x\n", x);
			return -1;
		}
	}

	// Print out the contents of the table (16 elements in each row)
	printf("forward\n");
	for (i = 0; i < 16; i++)
	{
		for (j = 0; j < 16; j++)
		{
			if (j == 15) printf("%x\n", forward((16 * i) + j));
			else printf("%x & ", forward((16 * i) + j));
		}
	}

	printf("inverse\n");
	for (i = 0; i < 16; i++)
	{
		for (j = 0; j < 16; j++)
		{
			if (j == 15) printf("%x\n", inverse((16 * i) + j));
			else printf("%x & ", inverse((16 * i) + j));
		}
	}

	// printf ("char S[256] = {\n");
	// for (i = 0; i < 16; i++) {
	// 	for (j = 0; j < 16; j++) {
	// 		printf ( "%3d, ", S[i*16+j]);
	// 	}
	// 	printf ( "\n" );
	// }
	// printf ( "};\n\n" );
	// printf ("char Si[256] = {\n");
	// for (i = 0; i < 16; i++) {
	// 	for (j = 0; j < 16; j++) {
	// 		printf ( "%3d, ", SINV[i*16+j]);
	// 	}
	// 	printf ( "\n" );
	// }
	// printf ( "};\n\n" );
}
