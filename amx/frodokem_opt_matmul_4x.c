#include <stdint.h>
#include <string.h>

#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifdef ALLOC_MMAP
uint16_t *tmp_opt_sa_4x;

__attribute__((constructor)) static void alloc_opt_sa_4x(void) {
    tmp_opt_sa_4x = MEMORY_ALLOC(32 * sizeof(uint16_t));
}
#endif

void amx_mul_add_sa_plus_e_up_to_32_cols_4x(uint16_t *out[4], const int16_t *A, const uint16_t *s, const uint16_t *e[4],
                                            int n, int col) {
    // The code assumes nbar == 8 throughout
    // Inputs: s', e' (N_BAR x N), A (N x N)
    // Output: out = s'*A + e' (N_BAR x N)

    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Z_NO_SKIP |
                                MAC16_Y_NO_SKIP | MAC16_X_NO_SKIP | MAC16_X_ENABLE_MODE(0) | MAC16_X_ENABLE_VALUE(0);

#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_opt_sa_4x;
#else
    uint16_t tmp[32];
#endif

    for (int b = 0; b < 4; b++) {
        for (int k = 0; k < 7; k++) {
            AMX_LDZ(AMX_PTR(&e[b][k * n + col]) | LDZ_Z_ROW(2 * (8 * b + k)) | LDZ_LOAD_SINGLE);
        }

        if (col + 32 <= n) {
            AMX_LDZ(AMX_PTR(&e[b][7 * n + col]) | LDZ_Z_ROW(2 * (8 * b + 7)) | LDZ_LOAD_SINGLE);
        }
        else {
            int sz = n - col;
            AMX_LDX(AMX_PTR(&e[b][8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - sz)) | MAC16_Y_SKIP | MAC16_Z_SKIP |
                      MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(sz) | MAC16_Z_ROW(2 * (8 * b + 7)));
        }
    }

    for (int j0 = 0; j0 < n; j0 += 32) {
        int stride = MIN((n - j0) / 4, 8);

        for (int j = j0; j < j0 + stride; j++) {
            AMX_LDY(AMX_PTR(&s[j * 32]) | LDY_REG(0) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(j + stride) * 32]) | LDY_REG(1) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(j + 2 * stride) * 32]) | LDY_REG(2) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(j + 3 * stride) * 32]) | LDY_REG(3) | LDY_LOAD_SINGLE);

            AMX_LDX(AMX_PTR(&A[j * 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + stride) * 32]) | LDX_REG(1) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 2 * stride) * 32]) | LDX_REG(2) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 3 * stride) * 32]) | LDX_REG(3) | LDX_LOAD_SINGLE);

            AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_REG(0));
            AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_REG(1));
            AMX_MAC16(mac16_args | MAC16_X_REG(2) | MAC16_Y_REG(2));
            AMX_MAC16(mac16_args | MAC16_X_REG(3) | MAC16_Y_REG(3));
        }
    }

    if (col + 32 <= n) {
        for (int b = 0; b < 4; b++) {
            for (int k = 0; k < 8; k++) {
                AMX_STZ(AMX_PTR(&out[b][k * n + col]) | STZ_Z_ROW(2 * (8 * b + k)) | STZ_STORE_SINGLE);
            }
        }
    }
    else {
        for (int b = 0; b < 4; b++) {
            for (int k = 0; k < 8; k++) {
                AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(2 * (8 * b + k)) | STZ_STORE_SINGLE);
                memcpy(&out[b][k * n + col], tmp, (n - col) * sizeof(uint16_t));
            }
        }
    }
}

void amx_mul_add_sa_plus_e_up_to_32_rows_4x(uint16_t *out[4], const int16_t *A, const uint16_t *s, const uint16_t *e[4],
                                            int n, int row) {
    // The code assumes nbar == 8 throughout
    // Inputs: s', e' (N_BAR x N), A (N x N)
    // Output: out = s'*A + e' (N_BAR x N)

    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Z_NO_SKIP |
                                MAC16_Y_NO_SKIP | MAC16_X_NO_SKIP | MAC16_X_ENABLE_MODE(0) | MAC16_X_ENABLE_VALUE(0);
    const int n32 = 32 * ((n + 31) / 32);

#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_opt_sa_4x;
#else
    uint16_t tmp[32];
#endif

    const uint16_t **in = (row == 0) ? e : (const uint16_t **)out;

    for (int i0 = 0; i0 < n; i0 += 32) {
        int stride = MIN((n - row) / 4, 8);

        for (int b = 0; b < 4; b++) {
            for (int k = 0; k < 7; k++) {
                AMX_LDZ(AMX_PTR(&in[b][k * n + i0]) | LDZ_Z_ROW(2 * (8 * b + k)) | LDZ_LOAD_SINGLE);
            }

            if (i0 + 32 <= n) {
                AMX_LDZ(AMX_PTR(&in[b][7 * n + i0]) | LDZ_Z_ROW(2 * (8 * b + 7)) | LDZ_LOAD_SINGLE);
            }
            else {
                int sz = n - i0;
                AMX_LDX(AMX_PTR(&in[b][8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
                AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - sz)) | MAC16_Y_SKIP | MAC16_Z_SKIP |
                          MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(sz) | MAC16_Z_ROW(2 * (8 * b + 7)));
            }
        }

        for (int j = 0; j < stride; j++) {
            AMX_LDY(AMX_PTR(&s[(row + j) * 32]) | LDY_REG(0) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(row + j + stride) * 32]) | LDY_REG(1) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(row + j + 2 * stride) * 32]) | LDY_REG(2) | LDY_LOAD_SINGLE);
            AMX_LDY(AMX_PTR(&s[(row + j + 3 * stride) * 32]) | LDY_REG(3) | LDY_LOAD_SINGLE);

            AMX_LDX(AMX_PTR(&A[j * n32 + i0]) | LDX_REG(0) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + stride) * n32 + i0]) | LDX_REG(1) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 2 * stride) * n32 + i0]) | LDX_REG(2) | LDX_LOAD_SINGLE);
            AMX_LDX(AMX_PTR(&A[(j + 3 * stride) * n32 + i0]) | LDX_REG(3) | LDX_LOAD_SINGLE);

            AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_REG(0));
            AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_REG(1));
            AMX_MAC16(mac16_args | MAC16_X_REG(2) | MAC16_Y_REG(2));
            AMX_MAC16(mac16_args | MAC16_X_REG(3) | MAC16_Y_REG(3));
        }

        if (i0 + 32 <= n) {
            for (int b = 0; b < 4; b++) {
                for (int k = 0; k < 8; k++) {
                    AMX_STZ(AMX_PTR(&out[b][k * n + i0]) | STZ_Z_ROW(2 * (8 * b + k)) | STZ_STORE_SINGLE);
                }
            }
        }
        else {
            for (int b = 0; b < 4; b++) {
                for (int k = 0; k < 8; k++) {
                    AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(2 * (8 * b + k)) | STZ_STORE_SINGLE);
                    memcpy(&out[b][k * n + i0], tmp, (n - i0) * sizeof(uint16_t));
                }
            }
        }
    }
}
