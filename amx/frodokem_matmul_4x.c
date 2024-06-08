#include <stdint.h>
#include <string.h>

#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

void amx_mul_add_sb_plus_e_4x(uint16_t *out[4], const uint16_t *b, const uint16_t *s[4], const uint16_t *e[4], int n) {
    // The code assumes nbar == 8 throughout
    // Multiply by s on the left
    // Inputs: b (N x N_BAR), s (N_BAR x N), e (N_BAR x N_BAR)
    // Output: out = s*b + e (N_BAR x N_BAR)
    const uint64_t mac16_args = MAC16_MATRIX | MAC16_Y_REG(0) | MAC16_X_REG(0) | MAC16_Z_ROW(0) | MAC16_Z_NO_SKIP |
                                MAC16_Y_NO_SKIP | MAC16_X_NO_SKIP | MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(8);
    const uint64_t extrv_args = EXTRV_GENERIC | EXTRV_DESTINATION_Y | EXTRV_LANE_WIDTH_16_BIT;

    for (int i = 0; i < 4; i++) {
        AMX_LDX(AMX_PTR(&e[i][0]) | LDX_REG(0) | LDX_LOAD_SINGLE);
        AMX_LDX(AMX_PTR(&e[i][32]) | LDX_REG(1) | LDX_LOAD_SINGLE);

        for (int k = 0; k < 8; k++) {
            AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * 8 * k) | MAC16_Y_SKIP | MAC16_Z_SKIP | MAC16_X_ENABLE_MODE(2) |
                      MAC16_X_ENABLE_VALUE(8) | MAC16_Z_ROW(2 * (8 * i + k)));
        }
    }

    for (int j0 = 0; j0 < n; j0 += 32) {
        if (j0 + 32 <= n) {
            for (int i = 0; i < 4; i++) {
                for (int k = 0; k < 8; k++) {
                    AMX_LDZ(AMX_PTR(&s[i][k * n + j0]) | LDZ_Z_ROW(2 * (8 * i + k) + 1) | LDZ_LOAD_SINGLE);
                }
            }

            AMX_LDX_QUAD(&b[j0 * 8], 0);
            AMX_LDX_QUAD(&b[(j0 + 16) * 8], 4);

            for (int j = j0; j < j0 + 8; j++) {
#ifdef CPU_M1
                AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * (j - j0) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (j - j0 + 8) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (j - j0 + 2 * 8) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (j - j0 + 3 * 8) + 1));
#else
                AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * (j - j0) + 1) | EXTRV_MULTIPLE_4);
#endif

                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0)) | MAC16_Y_REG(0));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 8)) | MAC16_Y_REG(1));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 2 * 8)) | MAC16_Y_REG(2));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 3 * 8)) | MAC16_Y_REG(3));
            }
        }
        else {
            for (int i = 0; i < 4; i++) {
                for (int k = 0; k < 7; k++) {
                    AMX_LDZ(AMX_PTR(&s[i][k * n + j0]) | LDZ_Z_ROW(2 * (8 * i + k) + 1) | LDZ_LOAD_SINGLE);
                }

                int sz = n - j0;
                AMX_LDX(AMX_PTR(&s[i][8 * n - 32]) | LDX_REG(0) | LDX_LOAD_SINGLE);
                AMX_MAC16(MAC16_VECTOR | MAC16_X_OFFSET(2 * (32 - sz)) | MAC16_Y_SKIP | MAC16_Z_SKIP |
                          MAC16_X_ENABLE_MODE(2) | MAC16_X_ENABLE_VALUE(sz) | MAC16_Z_ROW(2 * (8 * i + 7) + 1));
            }

            AMX_LDX_QUAD(&b[j0 * 8], 0);

            for (int j = j0; j < j0 + 4; j++) {
                AMX_EXTRV(extrv_args | EXTRV_REG(0) | EXTRV_Z_COLUMN(2 * (j - j0) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(1) | EXTRV_Z_COLUMN(2 * (j - j0 + 4) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(2) | EXTRV_Z_COLUMN(2 * (j - j0 + 2 * 4) + 1));
                AMX_EXTRV(extrv_args | EXTRV_REG(3) | EXTRV_Z_COLUMN(2 * (j - j0 + 3 * 4) + 1));

                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0)) | MAC16_Y_REG(0));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 4)) | MAC16_Y_REG(1));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 2 * 4)) | MAC16_Y_REG(2));
                AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * 8 * (j - j0 + 3 * 4)) | MAC16_Y_REG(3));
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 8; k++) {
            AMX_EXTRH(EXTRH_COPY_ONLY | EXTRH_COPY_ONLY_OFFSET(2 * 8 * k) | EXTRH_COPY_ONLY_Z_ROW(2 * (8 * i + k)) |
                      EXTRH_COPY_ONLY_LANE_WIDTH_16_BIT | EXTRH_COPY_ONLY_ENABLE_MODE(2) |
                      EXTRH_COPY_ONLY_ENABLE_VALUE(8));
        }

        AMX_STX(AMX_PTR(&out[i][0]) | STX_REG(0) | STX_STORE_SINGLE);
        AMX_STX(AMX_PTR(&out[i][32]) | STX_REG(1) | STX_STORE_SINGLE);
    }
}
