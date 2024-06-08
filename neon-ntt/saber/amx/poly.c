#include <stdint.h>
#include <string.h>

#include "SABER_params.h"
#include "aarch64.h"
#include "amx.h"
#include "polymodmul.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

// We use this function to initialize AMX. The use of __attribute__((constructor)) ensures this function is called
// before main().
__attribute__((constructor)) void init_AMX(void) {
    AMX_SET();
}

void VectorAdd(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N]) {
    const uint64_t vecint_args = VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5) | VECINT_ALU_MODE_Z_ADD_X_ADD_Y;

    AMX_LDX_QUAD(AMX_PTR(a), 0);
    AMX_LDX_QUAD(AMX_PTR(&a[128]), 4);

    AMX_LDZ(AMX_PTR(b) | LDZ_Z_ROW(0) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[32]) | LDZ_Z_ROW(16) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[64]) | LDZ_Z_ROW(32) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[96]) | LDZ_Z_ROW(48) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[128]) | LDZ_Z_ROW(8) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[160]) | LDZ_Z_ROW(24) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[192]) | LDZ_Z_ROW(40) | LDZ_LOAD_SINGLE);
    AMX_LDZ(AMX_PTR(&b[224]) | LDZ_Z_ROW(56) | LDZ_LOAD_SINGLE);

    for (int i = 0; i < 8; i += 4) {
#ifdef CPU_M1
        AMX_VECINT(vecint_args | VECINT_X_REG(i) | VECINT_Z_ROW(2 * i));
        AMX_VECINT(vecint_args | VECINT_X_REG(i + 1) | VECINT_Z_ROW(2 * i + 16));
        AMX_VECINT(vecint_args | VECINT_X_REG(i + 2) | VECINT_Z_ROW(2 * i + 32));
        AMX_VECINT(vecint_args | VECINT_X_REG(i + 3) | VECINT_Z_ROW(2 * i + 48));
#else
        AMX_VECINT(vecint_args | VECINT_X_REG(i) | VECINT_Z_ROW(2 * i) | VECINT_MULTIPLE_4);
#endif

        AMX_STZ(AMX_PTR(&res[32 * i]) | STZ_Z_ROW(2 * i) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[32 * (i + 1)]) | STZ_Z_ROW(2 * i + 16) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[32 * (i + 2)]) | STZ_Z_ROW(2 * i + 32) | STZ_STORE_SINGLE);
        AMX_STZ(AMX_PTR(&res[32 * (i + 3)]) | STZ_Z_ROW(2 * i + 48) | STZ_STORE_SINGLE);
    }
}

#ifdef ALLOC_MMAP
static uint16_t *tmp_ip;

__attribute__((constructor)) static void alloc_ip(void) {
    tmp_ip = MEMORY_ALLOC(SABER_N * sizeof(uint16_t));
}
#endif

void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]) {
#ifdef ALLOC_STACK
    uint16_t tmp[SABER_N];
#else
    uint16_t *tmp = tmp_ip;
#endif

    amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(res, b[0], s[0], 256, 8);

    for (int i = 1; i < SABER_L; i++) {
        amx_poly_mul_mod_65536_mod_x_d_plus_1_u16_32nx32n_coeffs(tmp, b[i], s[i], 256, 8);
        VectorAdd(res, tmp, res);
    }
}
