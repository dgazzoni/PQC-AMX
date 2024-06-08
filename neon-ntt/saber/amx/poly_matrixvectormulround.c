#include <stdint.h>
#include <string.h>

#include "SABER_params.h"
#include "aarch64.h"
#include "amx.h"
#include "poly.h"
#include "polymodmul.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#define h1 (1 << (SABER_EQ - SABER_EP - 1))

#ifdef ALLOC_MMAP
static uint16_t *tmp_mvm, *h1v_mvm;

__attribute__((constructor)) static void alloc_mvm(void) {
    tmp_mvm = MEMORY_ALLOC(SABER_N * sizeof(uint16_t));
    h1v_mvm = MEMORY_ALLOC(32 * sizeof(uint16_t));

    for (int i = 0; i < 32; i++) {
        h1v_mvm[i] = h1;
    }
}
#endif

void AddHShift(uint16_t res[SABER_L][SABER_N]) {
#ifdef ALLOC_MMAP
    uint16_t *h1v = h1v_mvm;
#else
    const uint16_t h1v[32] = {h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1,
                              h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1, h1};
#endif
    const uint64_t vecint_add_args = VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5) | VECINT_ALU_MODE_Z_ADD_X_ADD_Y;
    const uint64_t vecint_shift_args = VECINT_ALU_MODE_Z_SHIFT_RIGHT_S | VECINT_RIGHT_SHIFT_AMOUNT(SABER_EQ - SABER_EP);

    AMX_LDX(AMX_PTR(h1v) | LDX_REG(0) | LDX_LOAD_SINGLE);

    for (int i = 0; i < SABER_L; i++) {
        AMX_LDZ(AMX_PTR(&res[i][0]) | LDZ_Z_ROW(0) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][32]) | LDZ_Z_ROW(8) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][64]) | LDZ_Z_ROW(16) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][96]) | LDZ_Z_ROW(24) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][128]) | LDZ_Z_ROW(32) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][160]) | LDZ_Z_ROW(40) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][192]) | LDZ_Z_ROW(48) | LDZ_LOAD_SINGLE);
        AMX_LDZ(AMX_PTR(&res[i][224]) | LDZ_Z_ROW(56) | LDZ_LOAD_SINGLE);

        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(0));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(8));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(16));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(24));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(32));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(40));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(48));
        AMX_VECINT(vecint_add_args | VECINT_X_REG(0) | VECINT_Z_ROW(56));

#ifdef CPU_M1
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(0));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(8));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(16));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(24));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(32));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(40));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(48));
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(56));
#else
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(0) | VECINT_MULTIPLE_4);
        AMX_VECINT(vecint_shift_args | VECINT_Z_ROW(8) | VECINT_MULTIPLE_4);
#endif

        AMX_STZ(AMX_PTR(&res[i][0]) | STZ_Z_ROW(0) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][32]) | STZ_Z_ROW(8) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][64]) | STZ_Z_ROW(16) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][96]) | STZ_Z_ROW(24) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][128]) | STZ_Z_ROW(32) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][160]) | STZ_Z_ROW(40) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][192]) | STZ_Z_ROW(48) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[i][224]) | STZ_Z_ROW(56) | STZ_STORE_SINGLE);
    }
}

void MatrixVectorMulRound(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N],
                          uint16_t res[SABER_L][SABER_N], int16_t transpose) {
#ifdef ALLOC_STACK
    uint16_t tmp[SABER_N];
#else
    uint16_t *tmp = tmp_mvm;
#endif
    for (int i = 0; i < SABER_L; i++) {
        if (transpose == 1) {
            amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(res[i], A[0][i], s[0], 256, 8);
            for (int j = 1; j < SABER_L; j++) {
                amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(tmp, A[j][i], s[j], 256, 8);
                VectorAdd(res[i], tmp, res[i]);
            }
        }
        else {
            amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(res[i], A[i][0], s[0], 256, 8);
            for (int j = 1; j < SABER_L; j++) {
                amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(tmp, A[i][j], s[j], 256, 8);
                VectorAdd(res[i], tmp, res[i]);
            }
        }
    }

    AddHShift(res);
}
