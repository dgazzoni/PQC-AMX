#ifndef FRODOKEM_MATMUL_4X_H
#define FRODOKEM_MATMUL_4X_H

#include <stdint.h>

void amx_mul_add_sb_plus_e_4x(uint16_t *out[4], const uint16_t *b, const uint16_t *s[4], const uint16_t *e[4],
                              int n);

#endif  // FRODOKEM_MATMUL_4X_H
