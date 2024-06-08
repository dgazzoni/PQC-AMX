#ifndef POLY_H
#define POLY_H

#include <stdint.h>

#include "SABER_params.h"

void VectorAdd(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N]);
void MatrixVectorMulRound(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N],
                          uint16_t res[SABER_L][SABER_N], int16_t transpose);
void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);

#endif  // POLY_H
