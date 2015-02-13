/*
 * File: galois.h
 */

#ifndef GALOIS_H_
#define GALOIS_H_

#include <stdint.h>
#include <stdio.h>

// Quotient and remainder struct
typedef struct 
{
	uint16_t q;
	uint16_t r;
	uint8_t error;
} QR;

// Bit masks for the MSB and LSB
#define MSB_16  0x8000
#define HMSB_16 0x10000
#define MSB_8   0x80
#define HMSB_8  0x100
#define LSB     0x1

// Polynomial arithmetic in GF(2^8)
uint8_t g8_add(uint8_t x, uint8_t y);
uint8_t g8_mul(uint8_t x, uint8_t y);
QR g8_div(uint16_t ai, uint8_t b);
uint8_t g8_inv(uint8_t x);

// Polynomial arithmetic in GF(2^16)
uint16_t g16_add(uint16_t x, uint16_t y);
uint16_t g16_sub(uint16_t x, uint16_t y);
uint16_t g16_mul(uint16_t x, uint16_t y);
QR g16_div(uint32_t ai, uint16_t b);
uint16_t g16_inv(uint16_t x);

// Basis change functions (bit-matrix multiplication).
uint16_t g16_change_basis(uint16_t x, uint16_t* M);
uint8_t g8_change_basis(uint8_t x, uint8_t* M);

// Arithmetic in GF(2^2)
uint8_t g22_mul(uint8_t x, uint8_t y);
uint8_t g22_sq(uint8_t x); // same as inverse according to FLT

// Arithmetic in GF((2^2)^2)
uint8_t g222_mul(uint8_t x, uint8_t y);
uint8_t g222_sq(uint8_t x);
uint8_t g222_inv(uint8_t x);

// Arithmetic in GF(((2^2)^2)^2)
uint8_t g2222_inv(uint8_t x);
uint8_t g2222_mul(uint8_t x, uint8_t y);

// Arithmetic in GF((((2^2)^2)^2)^2)
uint16_t g22222_inv(uint16_t x);

#endif /* GALOIS_H_ */