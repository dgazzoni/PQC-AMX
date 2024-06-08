#ifndef FRODOKEM_OPT_MATMUL_4X_H
#define FRODOKEM_OPT_MATMUL_4X_H

#include <stdint.h>

void amx_mul_add_sa_plus_e_up_to_32_cols_4x(uint16_t *out[4], const int16_t *A, const uint16_t *s, const uint16_t *e[4],
                                            int n, int col);
void amx_mul_add_sa_plus_e_up_to_32_rows_4x(uint16_t *out[4], const int16_t *A, const uint16_t *s, const uint16_t *e[4],
                                            int n, int row);

#endif  // FRODOKEM_OPT_MATMUL_4X_H
