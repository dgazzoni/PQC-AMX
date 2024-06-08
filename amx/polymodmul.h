#ifndef POLYMODMUL_H
#define POLYMODMUL_H

#include <stdint.h>

void amx_poly_mul_mod_65536_mod_x_d_minus_1_u16_32nx32n_coeffs(uint16_t xy[], const uint16_t x[], const uint16_t y[],
                                                               int d, int n);
void amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(uint16_t xy[], const uint16_t x[], const uint16_t y[],
                                                              int d, int n);

#endif
