#include <stdint.h>
#include <string.h>

#include "SABER_params.h"
#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#ifdef ALLOC_MMAP
static uint16_t *tmp_PAXPY;
static int16_t *minus_one_PAXPY;

__attribute__((constructor)) void alloc_PAXPY(void) {
    tmp_PAXPY = MEMORY_ALLOC(32 * sizeof(uint16_t));
    minus_one_PAXPY = MEMORY_ALLOC(32 * sizeof(uint16_t));

    for (int i = 0; i < 32; i++) {
        minus_one_PAXPY[i] = -1;
    }
}
#else
static const int16_t minus_one_PAXPY[32] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
#endif

// This is algorithm III.1 in the paper, with some loop unrolling.
static void prepare_matrix_A(const uint16_t A[SABER_L][SABER_L][SABER_N], int col, int16_t transpose) {
    const uint64_t mac16_args = MAC16_VECTOR | MAC16_Z_SKIP;
    AMX_LDY(AMX_PTR(minus_one_PAXPY) | LDY_REG(0) | LDY_LOAD_SINGLE);

    for (int j = 0; j < SABER_L; j++) {
#pragma GCC unroll 2
        for (int i = 0; i < 8; i += 4) {
            int r, c;
            if (transpose == 1) {
                r = col;
                c = j;
            }
            else {
                r = j;
                c = col;
            }

            // Positive values: bottom rows

            AMX_LDZ(AMX_PTR(&A[r][c][32 * i]) | LDZ_Z_ROW(2 * (SABER_L * i + j + 8 * SABER_L) + 1) | LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 1)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 1) + j + 8 * SABER_L) + 1) |
                    LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 2)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 2) + j + 8 * SABER_L) + 1) |
                    LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 3)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 3) + j + 8 * SABER_L) + 1) |
                    LDZ_LOAD_SINGLE);

            AMX_LDX_QUAD(&A[r][c][32 * i], i);

            // Negative values: top rows

            if (i != 0) {
                AMX_MAC16(mac16_args | MAC16_X_REG(i) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (SABER_L * i + j) + 1));
            }
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 1) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (SABER_L * (i + 1) + j) + 1));
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 2) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (SABER_L * (i + 2) + j) + 1));
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 3) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (SABER_L * (i + 3) + j) + 1));
        }
    }
}

// This corresponds to the initial part of Algorithm III.2 in the paper.
static void prepare_matrix_B(const uint16_t sp[SABER_N], uint16_t tmp[32]) {
    AMX_LDX_QUAD(AMX_PTR(sp), 0);
    AMX_LDX_QUAD(AMX_PTR(&sp[128]), 4);

    // Recall that we didn't write to the first SABER_L odd rows, so we don't need to save them here
    AMX_MAC16(MAC16_VECTOR | MAC16_X_REG(7) | MAC16_Y_REG(0) | MAC16_Z_ROW(1) | MAC16_Z_SKIP);
    AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(1) | STZ_STORE_SINGLE);
}

// This corresponds to the final part of Algorithm III.2 in the paper, with some loop unrolling performed.
static void saber_TMVP(const uint16_t sp[SABER_N], uint16_t tmp[32]) {
    const uint64_t mac16_args =
        MAC16_MATRIX | MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE((SABER_L * 8) & 31) | MAC16_Z_ROW(0);
    const uint64_t extrv_args = EXTRV_GENERIC | EXTRV_DESTINATION_Y | EXTRV_LANE_WIDTH_16_BIT;

    for (int j = 0; j < 32; j += 8) {
        AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * j + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (j + 1) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (j + 2) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (j + 3) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(4) | EXTRV_Z_COLUMN(2 * (j + 4) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(5) | EXTRV_Z_COLUMN(2 * (j + 5) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(6) | EXTRV_Z_COLUMN(2 * (j + 6) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(7) | EXTRV_Z_COLUMN(2 * (j + 7) + 1));

        AMX_LDX(AMX_PTR(tmp) | LDX_REG(7));

        AMX_MAC16(mac16_args | MAC16_X_OFFSET((2 * (256 - j)) & 511) | MAC16_Y_OFFSET(2 * (0 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 1))) | MAC16_Y_OFFSET(2 * (32 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 2))) | MAC16_Y_OFFSET(2 * (64 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 3))) | MAC16_Y_OFFSET(2 * (96 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 4))) | MAC16_Y_OFFSET(2 * (128 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 5))) | MAC16_Y_OFFSET(2 * (160 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 6))) | MAC16_Y_OFFSET(2 * (192 + SABER_L * 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 7))) | MAC16_Y_OFFSET(2 * (224 + SABER_L * 8)));

        AMX_LDX(AMX_PTR(&sp[32 * 7]) | LDX_REG(7));

#pragma GCC unroll 7
        for (int i = 1; i < 8; i++) {
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - j)) | MAC16_Y_OFFSET(2 * (0 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 1))) |
                      MAC16_Y_OFFSET(2 * (32 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 2))) |
                      MAC16_Y_OFFSET(2 * (64 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 3))) |
                      MAC16_Y_OFFSET(2 * (96 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 4))) |
                      MAC16_Y_OFFSET(2 * (128 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 5))) |
                      MAC16_Y_OFFSET(2 * (160 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 6))) |
                      MAC16_Y_OFFSET(2 * (192 + SABER_L * (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 7))) |
                      MAC16_Y_OFFSET(2 * (224 + SABER_L * (8 - i))));
        }
    }
}

void amx_PAXPY(const uint16_t s[SABER_N], const uint16_t A[SABER_L][SABER_L][SABER_N], int col, int16_t transpose) {
#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_PAXPY;
#else
    uint16_t tmp[32];
#endif

    prepare_matrix_A(A, col, transpose);
    prepare_matrix_B(s, tmp);
    saber_TMVP(s, tmp);
}

#define h1 (1 << (SABER_EQ - SABER_EP - 1))

#ifdef ALLOC_MMAP
static uint16_t *h1v_mvm;

__attribute__((constructor)) static void alloc_mvm(void) {
    h1v_mvm = MEMORY_ALLOC(32 * sizeof(uint16_t));

    for (int i = 0; i < 32; i++) {
        h1v_mvm[i] = h1;
    }
}
#endif

// This is Algorithm III.3 in the paper.
void MatrixVectorMulRound(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N],
                          uint16_t res[SABER_L][SABER_N], int16_t transpose) {
#ifdef ALLOC_MMAP
    uint16_t *h1v = h1v_mvm;
#else
    const uint16_t h1v[32] = {h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                              h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1};
#endif

    AMX_LDX(AMX_PTR(h1v) | LDX_REG(0) | LDX_LOAD_SINGLE);
    AMX_MAC16(MAC16_MATRIX | MAC16_X_REG(0) | MAC16_Y_SKIP | MAC16_Z_SKIP | MAC16_Z_ROW(0));

    for (int j = 0; j < SABER_L; j++) {
        amx_PAXPY(s[j], A, j, transpose);
    }

    for (int j = 0; j < SABER_L; j++) {
        for (int i = 0; i < 8; i++) {
            AMX_VECINT(VECINT_Z_ROW(2 * (SABER_L * i + j)) | VECINT_ALU_MODE_Z_SHIFT_RIGHT_S |
                       VECINT_RIGHT_SHIFT_AMOUNT(SABER_EQ - SABER_EP));
            AMX_STZ(AMX_PTR(&res[j][32 * i]) | STZ_Z_ROW(2 * (SABER_L * i + j)));
        }
    }
}
