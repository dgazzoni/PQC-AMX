#ifndef FRODOKEM_OPT_MATMUL_H
#define FRODOKEM_OPT_MATMUL_H

#include <stdint.h>

void amx_mul_add_as_plus_e_up_to_32_rows(uint16_t *out, const int16_t *A, const uint16_t *sT, const uint16_t *e, int n,
                                         int row);
void amx_mul_add_sa_plus_e_up_to_32_cols(uint16_t *out, const int16_t *A, const uint16_t *s, const uint16_t *e, int n,
                                         int col);
void amx_mul_add_sa_plus_e_up_to_32_rows(uint16_t *out, const int16_t *A, const uint16_t *s, const uint16_t *e, int n,
                                         int row);

#endif  // FRODOKEM_OPT_MATMUL_H
