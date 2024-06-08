#include <arm_neon.h>
#include <stdint.h>
#include <string.h>

#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#ifdef ALLOC_MMAP
static uint16_t *tmp_matmul;
static int16_t *minus_one;

__attribute__((constructor)) void alloc_matmul(void) {
    tmp_matmul = MEMORY_ALLOC(32 * sizeof(uint16_t));
    minus_one = MEMORY_ALLOC(32 * sizeof(uint16_t));

    for (int i = 0; i < 32; i++) {
        minus_one[i] = -1;
    }
}
#else
static const int16_t minus_one[32] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
#endif

static void prepare_v(const uint16_t y[]) {
    const uint64_t mac16_args = MAC16_VECTOR | MAC16_Z_SKIP;
    AMX_LDY(AMX_PTR(minus_one) | LDY_REG(0) | LDY_LOAD_SINGLE);

    for (int i = 0; i < 8; i += 4) {
        AMX_LDZ(AMX_PTR(&y[32 * i]) | LDZ_Z_ROW(2 * (i + 8) + 1) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&y[32 * (i + 1)]) | LDZ_Z_ROW(2 * (i + 1 + 8) + 1) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&y[32 * (i + 2)]) | LDZ_Z_ROW(2 * (i + 2 + 8) + 1) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&y[32 * (i + 3)]) | LDZ_Z_ROW(2 * (i + 3 + 8) + 1) | LDZ_LOAD_SINGLE);

        AMX_LDX_QUAD(&y[32 * i], i);

        AMX_MAC16(mac16_args | MAC16_X_REG(i) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * i + 1));
        AMX_MAC16(mac16_args | MAC16_X_REG(i + 1) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (i + 1) + 1));
        AMX_MAC16(mac16_args | MAC16_X_REG(i + 2) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (i + 2) + 1));
        AMX_MAC16(mac16_args | MAC16_X_REG(i + 3) | MAC16_Y_REG(0) | MAC16_Z_ROW(2 * (i + 3) + 1));
    }
}

static void prepare_M(const uint16_t x[], uint16_t tmp[32]) {
    AMX_LDX_QUAD(AMX_PTR(x), 0);
    AMX_LDX_QUAD(AMX_PTR(&x[128]), 4);

    AMX_STY(AMX_PTR(tmp) | STY_REG(0) | STY_STORE_SINGLE);

    AMX_LDY(AMX_PTR(minus_one) | LDY_REG(0) | LDY_LOAD_SINGLE);

    AMX_MAC16(MAC16_VECTOR | MAC16_X_REG(7) | MAC16_Y_REG(0) | MAC16_Z_SKIP);
    AMX_LDY(AMX_PTR(tmp) | LDY_REG(0) | LDY_LOAD_SINGLE);

    AMX_STZ(AMX_PTR(tmp) | STZ_Z_ROW(0) | STZ_STORE_SINGLE);
}

static void saber_matmul(const uint16_t x[], uint16_t tmp[32]) {
    const uint64_t mac16_args =
        MAC16_MATRIX | MAC16_MATRIX_Y_ENABLE_MODE(2) | MAC16_MATRIX_Y_ENABLE_VALUE(8) | MAC16_Z_ROW(0);
    const uint64_t extrv_args = EXTRV_GENERIC | EXTRV_DESTINATION_Y | EXTRV_LANE_WIDTH_16_BIT;
    uint64_t skip = MAC16_Z_SKIP;

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

        AMX_MAC16(mac16_args | MAC16_X_OFFSET((2 * (256 - j)) & 511) | MAC16_Y_OFFSET(2 * (0 + 8)) | skip);
        skip = 0ULL;
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 1))) | MAC16_Y_OFFSET(2 * (32 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 2))) | MAC16_Y_OFFSET(2 * (64 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 3))) | MAC16_Y_OFFSET(2 * (96 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 4))) | MAC16_Y_OFFSET(2 * (128 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 5))) | MAC16_Y_OFFSET(2 * (160 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 6))) | MAC16_Y_OFFSET(2 * (192 + 8)));
        AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (256 - (j + 7))) | MAC16_Y_OFFSET(2 * (224 + 8)));

        AMX_LDX(AMX_PTR(&x[32 * 7]) | LDX_REG(7));

#pragma GCC unroll 7
        for (int i = 1; i < 8; i++) {
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - j)) | MAC16_Y_OFFSET(2 * (0 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 1))) | MAC16_Y_OFFSET(2 * (32 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 2))) | MAC16_Y_OFFSET(2 * (64 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 3))) | MAC16_Y_OFFSET(2 * (96 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 4))) | MAC16_Y_OFFSET(2 * (128 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 5))) | MAC16_Y_OFFSET(2 * (160 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 6))) | MAC16_Y_OFFSET(2 * (192 + (8 - i))));
            AMX_MAC16(mac16_args | MAC16_X_OFFSET(2 * (32 * i - (j + 7))) | MAC16_Y_OFFSET(2 * (224 + (8 - i))));
        }
    }
}

void AMX_saber_polymodmul_matmul(const uint16_t x[], const uint16_t y[]) {
#ifdef ALLOC_MMAP
    uint16_t *tmp = tmp_matmul;
#else
    uint16_t tmp[32];
#endif

    prepare_v(y);

    prepare_M(x, tmp);

    saber_matmul(x, tmp);
}

void amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(uint16_t xy[], const uint16_t x[], const uint16_t y[],
                                                              int d, int n) {
    (void)d;
    (void)n;

    AMX_MAC16(MAC16_MATRIX | MAC16_X_SKIP | MAC16_Y_SKIP | MAC16_Z_SKIP | MAC16_Z_ROW(0));

    AMX_saber_polymodmul_matmul(x, y);

    for (int i = 0; i < 8; i++) {
        AMX_STZ(AMX_PTR(&xy[32 * i]) | STZ_Z_ROW(2 * i));
    }
}
