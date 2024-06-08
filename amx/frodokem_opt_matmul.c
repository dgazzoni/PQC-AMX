#include <stdint.h>
#include <string.h>

#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// We use this function to initialize AMX. The use of __attribute__((constructor)) ensures this function is called
// before main().
__attribute__((constructor)) void init_AMX(void) {
    AMX_SET();
}

#ifdef ALLOC_MMAP
uint16_t *tmp_opt_as;

__attribute__((constructor)) static void alloc_opt_as(void) {
    tmp_opt_as = MEMORY_ALLOC(32 * sizeof(uint16_t));
}
#endif

// This is algorithm IV.1 in the paper, with some modifications to handle the FrodoKEM-976 parameter set, in which
// n is not a multiple of 32; indeed, it is congruent to 16 (mod 32). We also perform some loop unrolling, so that we
// can use a special feature in the M2 and later SoCs, which is to execute certain operations multiple times using a
// single instruction.
void amx_mul_add_as_plus_e_up_to_32_rows(uint16_t *out, const int16_t *A, const uint16_t *sT, const uint16_t *e, int n,
                                         int row) {
    // The code assumes nbar == 8 throughout

    // Inputs: s, e (N x N_BAR), A (N x N)
    // Output: out = A*s + e (N x N_BAR)

    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Y_NO_SKIP |
                                MAC16_X_NO_SKIP | MAC16_Z_NO_SKIP | MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(8);
    const uint64_t extrv_args = EXTRV_GENERIC | EXTRV_DESTINATION_Y | EXTRV_LANE_WIDTH_16_BIT;
    const int n32 = 32 * ((n + 31) / 32);
#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_opt_as;
#else
    uint16_t tmp[32];
#endif

    for (int i = 0; i < MIN(n - row, 32); i++) {
        if ((row + i) < n - 3) {
            AMX_LDZ(AMX_PTR(&e[(row + i) * 8]) | LDZ_Z_ROW(2 * i) | LDZ_LOAD_SINGLE);
        }
        else {
            AMX_LDX(AMX_PTR(&e[8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - 8 * (n - (row + i)))) | MAC16_Y_SKIP | MAC16_Z_SKIP |
                      MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(8) | MAC16_Z_ROW(2 * i));
        }
    }

    for (int j0 = 0; j0 < n; j0 += 32) {
        for (int k = 0; k < MIN((n - j0) / 4, 8); k++) {
            AMX_LDX(AMX_PTR(&sT[j0 * 8 + k * 32]) | LDX_REG(k) | LDX_LOAD_SINGLE);
        }

        for (int i = 0; i < MIN(n - row, 32); i++) {
            AMX_LDZ(AMX_PTR(&A[i * n32 + j0]) | LDZ_Z_ROW(2 * i + 1) | LDZ_LOAD_SINGLE);
        }

        for (int j = j0; j < j0 + 8; j++) {
#ifdef CPU_M1
            AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * (j - j0) + 1));
            AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (j - j0 + 8) + 1));
            AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (j - j0 + 16) + 1));
            AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (j - j0 + 24) + 1));
#else
            AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * (j - j0) + 1) | EXTRV_MULTIPLE_4);
#endif

            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0)) | MAC16_Y_REG(0) | MAC16_MATRIX_Y_ENABLE_MODE(2) |
                      MAC16_MATRIX_Y_ENABLE_VALUE(MIN(n - row, 32) & 31));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 8)) | MAC16_Y_REG(1) |
                      MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE(MIN(n - row, 32) & 31));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 16)) | MAC16_Y_REG(2) |
                      MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE(MIN(n - row, 32) & 31));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 24)) | MAC16_Y_REG(3) |
                      MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE(MIN(n - row, 32) & 31));
        }
    }

    for (int i = 0; i < MIN(n - row, 32); i++) {
        if ((row + i) < n - 3) {
            AMX_STZ(AMX_PTR(&out[(row + i) * 8]) | STZ_Z_ROW(2 * i) | STZ_STORE_SINGLE);
        }
        else {
            AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(2 * i) | STZ_STORE_SINGLE);
            memcpy(&out[(row + i) * 8], tmp, 8 * sizeof(uint16_t));
        }
    }
}

#ifdef ALLOC_MMAP
uint16_t *tmp_opt_sa;

__attribute__((constructor)) static void alloc_opt_sa(void) {
    tmp_opt_sa = MEMORY_ALLOC(32 * sizeof(uint16_t));
}
#endif

void amx_mul_add_sa_plus_e_up_to_32_cols(uint16_t *out, const int16_t *A, const uint16_t *s, const uint16_t *e, int n,
                                         int col) {
    // The code assumes nbar == 8 throughout
    // Inputs: s', e' (N_BAR x N), A (N x N)
    // Output: out = s'*A + e' (N_BAR x N)

    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Z_NO_SKIP |
                                MAC16_Y_NO_SKIP | MAC16_X_NO_SKIP | MAC16_MATRIX_Y_ENABLE_MODE(2) |
                                MAC16_MATRIX_Y_ENABLE_VALUE(8) | MAC16_X_ENABLE_MODE(0) | MAC16_X_ENABLE_VALUE(0);

#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_opt_sa;
#else
    uint16_t tmp[32];
#endif

    for (int k = 0; k < 7; k++) {
        AMX_LDZ(AMX_PTR(&e[k * n + col]) | LDZ_Z_ROW(2 * k) | LDZ_LOAD_SINGLE);
    }

    if (col + 32 <= n) {
        AMX_LDZ(AMX_PTR(&e[7 * n + col]) | LDZ_Z_ROW(2 * 7) | LDZ_LOAD_SINGLE);
    }
    else {
        int sz = n - col;
        AMX_LDX(AMX_PTR(&e[8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
        AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - sz)) | MAC16_Y_SKIP | MAC16_Z_SKIP | MAC16_X_ENABLE_MODE(2) |
                  MAC16_X_ENABLE_VALUE(sz) | MAC16_Z_ROW(2 * 7));
    }

    for (int j0 = 0; j0 < n; j0 += 32) {
        int stride = MIN((n - j0) / 4, 8);

        for (int k = 0; k < stride; k++) {
            AMX_LDY(AMX_PTR(&s[j0 * 8 + k * 32]) | LDY_REG(k) | LDY_LOAD_SINGLE);
        }

        for (int j = j0; j < j0 + stride; j++) {
            AMX_LDX(AMX_PTR(&A[j * 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + stride) * 32]) | LDX_REG(1) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 2 * stride) * 32]) | LDX_REG(2) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 3 * stride) * 32]) | LDX_REG(3) | LDX_LOAD_SINGLE);

            AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_OFFSET(2 * 8 * (j - j0)));
            AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_OFFSET(2 * 8 * (j - j0 + stride)));
            AMX_MAC16(mac16_args | MAC16_X_REG(2) | MAC16_Y_OFFSET(2 * 8 * (j - j0 + 2 * stride)));
            AMX_MAC16(mac16_args | MAC16_X_REG(3) | MAC16_Y_OFFSET(2 * 8 * (j - j0 + 3 * stride)));
        }
    }

    if (col + 32 <= n) {
        for (int k = 0; k < 8; k++) {
            AMX_STZ(AMX_PTR(&out[k * n + col]) | STZ_Z_ROW(2 * k) | STZ_STORE_SINGLE);
        }
    }
    else {
        for (int k = 0; k < 8; k++) {
            AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(2 * k) | STZ_STORE_SINGLE);
            memcpy(&out[k * n + col], tmp, (n - col) * sizeof(uint16_t));
        }
    }
}

void amx_mul_add_sa_plus_e_up_to_32_rows(uint16_t *out, const int16_t *A, const uint16_t *s, const uint16_t *e, int n,
                                         int row) {
    // The code assumes nbar == 8 throughout
    // Inputs: s', e' (N_BAR x N), A (N x N)
    // Output: out = s'*A + e' (N_BAR x N)

    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Z_NO_SKIP |
                                MAC16_Y_NO_SKIP | MAC16_X_NO_SKIP | MAC16_MATRIX_Y_ENABLE_MODE(2) |
                                MAC16_MATRIX_Y_ENABLE_VALUE(8) | MAC16_X_ENABLE_MODE(0) | MAC16_X_ENABLE_VALUE(0);
    const int n32 = 32 * ((n + 31) / 32);

#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_opt_sa;
#else
    uint16_t tmp[32];
#endif

    const uint16_t *in = (row == 0) ? e : out;

    for (int i0 = 0; i0 < n; i0 += 32) {
        int stride = MIN((n - row) / 4, 8);

        for (int k = 0; k < 7; k++) {
            AMX_LDZ(AMX_PTR(&in[k * n + i0]) | LDZ_Z_ROW(2 * k) | LDZ_LOAD_SINGLE);
        }

        if (i0 + 32 <= n) {
            AMX_LDZ(AMX_PTR(&in[7 * n + i0]) | LDZ_Z_ROW(2 * 7) | LDZ_LOAD_SINGLE);
        }
        else {
            int sz = n - i0;
            AMX_LDX(AMX_PTR(&in[8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - sz)) | MAC16_Y_SKIP | MAC16_Z_SKIP |
                      MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(sz) | MAC16_Z_ROW(2 * 7));
        }

        for (int k = 0; k < stride; k++) {
            AMX_LDY(AMX_PTR(&s[row * 8 + k * 32]) | LDY_REG(k) | LDY_LOAD_SINGLE);
        }

        for (int j = 0; j < stride; j++) {
            AMX_LDX(AMX_PTR(&A[j * n32 + i0]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + stride) * n32 + i0]) | LDX_REG(1) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 2 * stride) * n32 + i0]) | LDX_REG(2) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 3 * stride) * n32 + i0]) | LDX_REG(3) | LDX_LOAD_SINGLE);

            AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_OFFSET(2 * 8 * j));
            AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_OFFSET(2 * 8 * (j + stride)));
            AMX_MAC16(mac16_args | MAC16_X_REG(2) | MAC16_Y_OFFSET(2 * 8 * (j + 2 * stride)));
            AMX_MAC16(mac16_args | MAC16_X_REG(3) | MAC16_Y_OFFSET(2 * 8 * (j + 3 * stride)));
        }

        if (i0 + 32 <= n) {
            for (int k = 0; k < 8; k++) {
                AMX_STZ(AMX_PTR(&out[k * n + i0]) | STZ_Z_ROW(2 * k) | STZ_STORE_SINGLE);
            }
        }
        else {
            for (int k = 0; k < 8; k++) {
                AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(2 * k) | STZ_STORE_SINGLE);
                memcpy(&out[k * n + i0], tmp, (n - i0) * sizeof(uint16_t));
            }
        }
    }
}
