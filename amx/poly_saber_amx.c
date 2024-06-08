#include <stdint.h>
#include <string.h>

#include "SABER_params.h"
#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#if SABER_L == 2
#define UNROLL_SABER_L _Pragma("GCC unroll 2")
#elif SABER_L == 3
#define UNROLL_SABER_L _Pragma("GCC unroll 3")
#elif SABER_L == 4
#define UNROLL_SABER_L _Pragma("GCC unroll 4")
#else
#error Invalid SABER_L, must be one of 2, 3, 4
#endif

#ifdef ALLOC_MMAP
static uint16_t *tmp_PAXPY, *v_PAXPY;
static int16_t *minus_one_PAXPY;

__attribute__((constructor)) void alloc_PAXPY(void) {
    tmp_PAXPY = MEMORY_ALLOC(32 * sizeof(uint16_t));
    v_PAXPY = MEMORY_ALLOC(2 * 4 * SABER_N * sizeof(uint16_t));
    minus_one_PAXPY = MEMORY_ALLOC(32 * sizeof(uint16_t));

    for (int i = 0; i < 32; i++) {
        minus_one_PAXPY[i] = -1;
    }
}
#else
static const int16_t minus_one_PAXPY[32] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
#endif

// This is algorithm III.1 in the paper, adapted for the Saber and FireSaber parameter sets. As discussed in the paper,
// unlike the LightSaber parameter set, the A matrix no longer fits in the odd-numbered Z registers. As such, they are
// preprocessed here and saved to the array "v", to be loaded on demand in the saber_mat_vec_mul() function. We also
// unroll some of the loops, especially to take advantage of certain instructions (such as EXTRV) which can be run
// multiple times using a single instruction, in the M2 and later SoCs.
static void prepare_matrix_A(const uint16_t A[SABER_L][SABER_L][SABER_N], uint16_t v[], int col, int16_t transpose) {
    const uint64_t mac16_args = MAC16_VECTOR | MAC16_Z_SKIP;
    const uint64_t extrv_args = EXTRV_GENERIC | EXTRV_DESTINATION_Y | EXTRV_LANE_WIDTH_16_BIT;

    UNROLL_SABER_L
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

            AMX_LDZ(AMX_PTR(&A[r][c][32 * i]) | LDZ_Z_ROW(2 * (SABER_L * i + j) + 1) | LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 1)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 1) + j) + 1) | LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 2)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 2) + j) + 1) | LDZ_LOAD_SINGLE);
            AMX_LDZ(AMX_PTR(&A[r][c][32 * (i + 3)]) | LDZ_Z_ROW(2 * (SABER_L * (i + 3) + j) + 1) | LDZ_LOAD_SINGLE);
        }
    }

    for (int i = 0; i < 8; i++) {
#ifdef CPU_M1
        AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * i + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (i + 8) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (i + 16) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (i + 24) + 1));
#else
        AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * i + 1) | EXTRV_MULTIPLE_4);
#endif

        AMX_STY(AMX_PTR(&v[(i + 32) * 32]) | STY_REG(0) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 8 + 32) * 32]) | STY_REG(1) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 16 + 32) * 32]) | STY_REG(2) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 24 + 32) * 32]) | STY_REG(3) | STZ_STORE_SINGLE);
    }

    AMX_LDY(AMX_PTR(minus_one_PAXPY) | LDY_REG(0) | LDY_LOAD_SINGLE);

    UNROLL_SABER_L
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

            AMX_LDX_QUAD(&A[r][c][32 * i], i);

            if (i != 0) {
                AMX_MAC16(mac16_args | MAC16_X_REG(i) | MAC16_Y_REG(0) |
                          MAC16_Z_ROW(2 * (8 * (4 - SABER_L) + SABER_L * i + j) + 1));
            }
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 1) | MAC16_Y_REG(0) |
                      MAC16_Z_ROW(2 * (8 * (4 - SABER_L) + SABER_L * (i + 1) + j) + 1));
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 2) | MAC16_Y_REG(0) |
                      MAC16_Z_ROW(2 * (8 * (4 - SABER_L) + SABER_L * (i + 2) + j) + 1));
            AMX_MAC16(mac16_args | MAC16_X_REG(i + 3) | MAC16_Y_REG(0) |
                      MAC16_Z_ROW(2 * (8 * (4 - SABER_L) + SABER_L * (i + 3) + j) + 1));
        }
    }

    for (int i = 0; i < 8; i++) {
#ifdef CPU_M1
        AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * i + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (i + 8) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (i + 16) + 1));
        AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (i + 24) + 1));
#else
        AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * i + 1) | EXTRV_MULTIPLE_4);
#endif
        AMX_STY(AMX_PTR(&v[(i)*32]) | STY_REG(0) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 8) * 32]) | STY_REG(1) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 16) * 32]) | STY_REG(2) | STZ_STORE_SINGLE);
        AMX_STY(AMX_PTR(&v[(i + 24) * 32]) | STY_REG(3) | STZ_STORE_SINGLE);
    }
}

// This corresponds to the initial part of Algorithm III.2 in the paper.
static void prepare_matrix_B(const uint16_t sp[SABER_N], uint16_t tmp[32]) {
    AMX_LDX_QUAD(AMX_PTR(sp), 0);
    AMX_LDX_QUAD(AMX_PTR(&sp[128]), 4);

    AMX_LDY(AMX_PTR(minus_one_PAXPY) | LDY_REG(0) | LDY_LOAD_SINGLE);

    AMX_MAC16(MAC16_VECTOR | MAC16_X_REG(7) | MAC16_Y_REG(0) | MAC16_Z_ROW(1) | MAC16_Z_SKIP);

    AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(1) | STZ_STORE_SINGLE);
}

// This corresponds to the final part of Algorithm III.2 in the paper, adapted for the Saber and FireSaber parameter
// sets. As the "A" matrix doesn't fit in the odd-numbered Z registers, we need to load the preprocessed values (from
// the prepare_matrix_A() function) on demand. We also perform some loop unrolling.
static void saber_TMVP(const uint16_t sp[SABER_N], uint16_t v[], uint16_t tmp[32]) {
    const uint64_t mac16_args =
        MAC16_MATRIX | MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE((SABER_L * 8) & 31) | MAC16_Z_ROW(0);

    for (int j = 0; j < 32; j += 4) {
        for (int k = j; k < j + 4; k++) {
            AMX_LDY(AMX_PTR(&v[k * 32]) | LDY_REG(2 * (k - j)) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&v[(k + 32) * 32]) | LDY_REG(2 * (k - j) + 1) | LDY_LOAD_SINGLE);
        }

        AMX_LDX(AMX_PTR(tmp) | LDX_REG(7));

        AMX_MAC16(mac16_args | MAC16_X_OFFSET((2 * (256 - j)) & 511) | MAC16_Y_REG(1));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 1))) | MAC16_Y_REG(3));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 2))) | MAC16_Y_REG(5));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 3))) | MAC16_Y_REG(7));

        AMX_LDX(AMX_PTR(&sp[32 * 7]) | LDX_REG(7));

#pragma GCC unroll 7
        for (int i = 1; i < 8; i++) {
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - j)) | MAC16_Y_OFFSET(2 * (0 + 32 - SABER_L * i)));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 1))) |
                      MAC16_Y_OFFSET(2 * (64 + 32 - SABER_L * i)));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 2))) |
                      MAC16_Y_OFFSET(2 * (128 + 32 - SABER_L * i)));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 3))) |
                      MAC16_Y_OFFSET(2 * (192 + 32 - SABER_L * i)));
        }
    }
}

void amx_PAXPY(const uint16_t s[SABER_N], const uint16_t A[SABER_L][SABER_L][SABER_N], int col, int16_t transpose) {
#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_PAXPY;
    uint16_t *v = v_PAXPY;
#else
    uint16_t tmp[32];
    uint16_t v[2 * 4 * SABER_N];
#endif

    prepare_matrix_A(A, v, col, transpose);
    prepare_matrix_B(s, tmp);
    saber_TMVP(s, v, tmp);
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
            AMX_VECINT(VECINT_X_REG(0) | VECINT_Z_ROW(2 * (SABER_L * i + j)) | VECINT_ALU_MODE_Z_SHIFT_RIGHT_S |
                       VECINT_RIGHT_SHIFT_AMOUNT(SABER_EQ - SABER_EP));
            AMX_STZ(AMX_PTR(&res[j][32 * i]) | STZ_Z_ROW(2 * (SABER_L * i + j)));
        }
    }
}
