#ifndef POLYMODMUL_MATMUL_H
#define POLYMODMUL_MATMUL_H

#include <arm_neon.h>
#include <stdint.h>

void negate_polynomial_slice_32(uint16_t res[], const uint16_t in[]);
void neon_transpose_8x8_block(uint16x8_t vr[]);
void prepare_toeplitz_matrix(uint16_t res[], const uint16_t in[]);
void AMX_matmul_8x32_by_32x32(uint16_t res[], const uint16_t BT[], const uint16_t A[], uint64_t skip_Z);
void AMX_saber_matmul_8x32_by_32x32(const uint16_t B_aux[], const uint16_t A_aux[], int index);
void AMX_saber_polymodmul_matmul(uint16_t res[], const uint16_t x[], const uint16_t y[]);

#endif  // POLYMODMUL_MATMUL_H
