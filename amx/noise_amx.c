#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "aarch64.h"
#include "amx.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

// From PQCrypto-LWEKE/FrodoKEM/src/frodo640.c, with some changes
#if PARAMS_N == 640
static const uint16_t AMX_CDF_TABLE_[32] = {0,
                                            1 + 4643,
                                            1 + 13363,
                                            1 + 20579,
                                            1 + 25843,
                                            1 + 29227,
                                            1 + 31145,
                                            1 + 32103,
                                            1 + 32525,
                                            1 + 32689,
                                            1 + 32745,
                                            1 + 32762,
                                            1 + 32766,
                                            1 + 32767,
                                            32768 + 1 + 4643,
                                            32768 + 1 + 13363,
                                            32768 + 1 + 20579,
                                            32768 + 1 + 25843,
                                            32768 + 1 + 29227,
                                            32768 + 1 + 31145,
                                            32768 + 1 + 32103,
                                            32768 + 1 + 32525,
                                            32768 + 1 + 32689,
                                            32768 + 1 + 32745,
                                            32768 + 1 + 32762,
                                            32768 + 1 + 32766};
static int16_t iota_[32] = {0,  1,  2,  3,  4,  5,  6,  7,   8,   9,   10,  11,  12,  0,   -1,  -2,
                            -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -12, -12, -12, -12, -12, -12};
#elif PARAMS_N == 976
static const uint16_t AMX_CDF_TABLE_[32] = {0,
                                            1 + 5638,
                                            1 + 15915,
                                            1 + 23689,
                                            1 + 28571,
                                            1 + 31116,
                                            1 + 32217,
                                            1 + 32613,
                                            1 + 32731,
                                            1 + 32760,
                                            1 + 32766,
                                            1 + 32767,
                                            32768 + 1 + 5638,
                                            32768 + 1 + 15915,
                                            32768 + 1 + 23689,
                                            32768 + 1 + 28571,
                                            32768 + 1 + 31116,
                                            32768 + 1 + 32217,
                                            32768 + 1 + 32613,
                                            32768 + 1 + 32731,
                                            32768 + 1 + 32760,
                                            32768 + 1 + 32766};
static int16_t iota_[32] = {0,  1,  2,  3,  4,  5,   6,   7,   8,   9,   10,  0,   -1,  -2,  -3,  -4,
                            -5, -6, -7, -8, -9, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10};
#elif PARAMS_N == 1344
static const uint16_t AMX_CDF_TABLE_[32] = {0,
                                            1 + 9142,
                                            1 + 23462,
                                            1 + 30338,
                                            1 + 32361,
                                            1 + 32725,
                                            1 + 32765,
                                            1 + 32767,
                                            32768 + 1 + 9142,
                                            32768 + 1 + 23462,
                                            32768 + 1 + 30338,
                                            32768 + 1 + 32361,
                                            32768 + 1 + 32725,
                                            32768 + 1 + 32765};
static int16_t iota_[32] = {0,  1,  2,  3,  4,  5,  6,  0,  -1, -2, -3, -4, -5, -6, -6, -6,
                            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6};
#else
#error Invalid/undefined PARAMS_N
#endif

static uint16_t shl_15_[32] = {32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768,
                               32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768,
                               32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768, 32768};

#ifdef ALLOC_MMAP
static uint16_t *AMX_CDF_TABLE, *shl_15;
static int16_t *iota;

__attribute__((constructor)) static void alloc_sample(void) {
    AMX_CDF_TABLE = (uint16_t *)MEMORY_ALLOC(32 * sizeof(uint16_t));
    iota = (int16_t *)MEMORY_ALLOC(32 * sizeof(uint16_t));
    shl_15 = (uint16_t *)MEMORY_ALLOC(32 * sizeof(uint16_t));

    memcpy(AMX_CDF_TABLE, AMX_CDF_TABLE_, 32 * sizeof(uint16_t));
    memcpy(iota, iota_, 32 * sizeof(uint16_t));
    memcpy(shl_15, shl_15_, 32 * sizeof(uint16_t));
}
#endif

// This is the loop body of Algorithm IV.2, unrolled to process 64 instead of 32 elements at a time. With the unrolling,
// it's possible to use a mode of certain instructions, such as "vecint", which allow processing multiple operations
// using a single instruction (in this case, 2 operations at a time), in the M2 and later SoCs.
void frodo_sample_n_64(uint16_t *s, const size_t n) {
#ifdef CPU_M1
    const uint64_t mac16_args = MAC16_VECTOR | MAC16_Z_SKIP;
    const uint64_t extrh_args = EXTRH_COPY_ONLY | EXTRH_COPY_ONLY_LANE_WIDTH_16_BIT;
#endif
    const uint64_t vecint_args = VECINT_ALU_MODE_Z_ADD_X_ADD_Y | VECINT_RIGHT_SHIFT_AMOUNT(1);
    const uint64_t genlut_generate_args = GENLUT_MODE_GENERATE_U16_U5 | GENLUT_SRC_X | GENLUT_DST_X | GENLUT_TABLE_Y;
    const uint64_t genlut_lookup_args = GENLUT_MODE_LOOKUP_16_U5 | GENLUT_SRC_X | GENLUT_DST_X | GENLUT_TABLE_Y;
    unsigned int i;
#ifdef ALLOC_STACK
    const uint16_t *AMX_CDF_TABLE = AMX_CDF_TABLE_, *shl_15 = shl_15_;
    const int16_t *iota = iota_;
#endif

    AMX_LDY(AMX_PTR(AMX_CDF_TABLE) | LDY_REG(0) | LDX_LOAD_SINGLE);
    AMX_LDY(AMX_PTR(iota) | LDX_REG(1) | LDX_LOAD_SINGLE);
    AMX_LDY(AMX_PTR(shl_15) | LDY_REG(2) | LDY_LOAD_SINGLE);

    for (i = 0; i < n; i += 64) {
        AMX_LDX(AMX_PTR(&s[i]) | LDX_REG(0) | LDX_LOAD_PAIR);

        // The next few instructions perform a 16-bit rotate right of s[i] by 1

#ifdef CPU_M1
        AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_REG(2) | MAC16_Z_ROW(0));
        AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_REG(2) | MAC16_Z_ROW(32));
#else
        AMX_VECINT(VECINT_X_REG(0) | VECINT_Y_REG(2) | VECINT_Z_ROW(0) | VECINT_ALU_MODE_X_MUL_Y |
                   VECINT_BROADCAST_MODE(3) | VECINT_MULTIPLE_2);
#endif

#ifdef CPU_M1
        AMX_VECINT(vecint_args | VECINT_X_REG(0) | VECINT_Z_ROW(0) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
        AMX_VECINT(vecint_args | VECINT_X_REG(1) | VECINT_Z_ROW(32) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
#else
        AMX_VECINT(vecint_args | VECINT_X_REG(0) | VECINT_Z_ROW(0) | VECINT_BROADCAST_MODE(5) | VECINT_MULTIPLE_2);
#endif

#ifdef CPU_M1
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(0) | EXTRH_COPY_ONLY_Z_ROW(0));
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(1) | EXTRH_COPY_ONLY_Z_ROW(32));
#else
        AMX_EXTRH(EXTRH_GENERIC | EXTRH_REG(0) | EXTRH_DESTINATION_X | EXTRH_Z_ROW(0) | EXTRH_LANE_WIDTH_16_BIT |
                  EXTRH_MULTIPLE_2);
#endif

        // Table lookup
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(0) | GENLUT_DST_REG(0) | GENLUT_TABLE_REG(0));
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(1) | GENLUT_DST_REG(1) | GENLUT_TABLE_REG(0));

        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(0) | GENLUT_DST_REG(0) | GENLUT_TABLE_REG(1));
        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(1) | GENLUT_DST_REG(1) | GENLUT_TABLE_REG(1));

        AMX_STX(AMX_PTR(&s[i]) | STX_REG(0) | STX_STORE_PAIR);
    }
}

// This is the loop body of Algorithm IV.2, unrolled to process 128 instead of 32 elements at a time. With the
// unrolling, it's possible to use a mode of certain instructions, such as "vecint", which allow processing multiple operations
// using a single instruction (in this case, 4 operations at a time), in the M2 and later SoCs.
void frodo_sample_n_128(uint16_t *s, const size_t n) {
#ifdef CPU_M1
    const uint64_t mac16_args = MAC16_VECTOR | MAC16_Z_SKIP;
    const uint64_t extrh_args = EXTRH_COPY_ONLY | EXTRH_COPY_ONLY_LANE_WIDTH_16_BIT;
#endif
    const uint64_t vecint_args = VECINT_ALU_MODE_Z_ADD_X_ADD_Y | VECINT_RIGHT_SHIFT_AMOUNT(1);
    const uint64_t genlut_generate_args = GENLUT_MODE_GENERATE_U16_U5 | GENLUT_SRC_X | GENLUT_DST_X | GENLUT_TABLE_Y;
    const uint64_t genlut_lookup_args = GENLUT_MODE_LOOKUP_16_U5 | GENLUT_SRC_X | GENLUT_DST_X | GENLUT_TABLE_Y;
    unsigned int i;
#ifdef ALLOC_STACK
    const uint16_t *AMX_CDF_TABLE = AMX_CDF_TABLE_, *shl_15 = shl_15_;
    const int16_t *iota = iota_;
#endif

    AMX_LDY(AMX_PTR(AMX_CDF_TABLE) | LDY_REG(0) | LDX_LOAD_SINGLE);
    AMX_LDY(AMX_PTR(iota) | LDY_REG(1) | LDY_LOAD_SINGLE);
    AMX_LDY(AMX_PTR(shl_15) | LDY_REG(2) | LDY_LOAD_SINGLE);

    for (i = 0; i < n; i += 128) {
        AMX_LDX_QUAD(AMX_PTR(&s[i]), 0);

        // The next few instructions perform a 16-bit rotate right of s[i] by 1

#ifdef CPU_M1
        AMX_MAC16(mac16_args | MAC16_X_REG(0) | MAC16_Y_REG(2) | MAC16_Z_ROW(0));
        AMX_MAC16(mac16_args | MAC16_X_REG(1) | MAC16_Y_REG(2) | MAC16_Z_ROW(16));
        AMX_MAC16(mac16_args | MAC16_X_REG(2) | MAC16_Y_REG(2) | MAC16_Z_ROW(32));
        AMX_MAC16(mac16_args | MAC16_X_REG(3) | MAC16_Y_REG(2) | MAC16_Z_ROW(48));
#else
        AMX_VECINT(VECINT_X_REG(0) | VECINT_Y_REG(2) | VECINT_Z_ROW(0) | VECINT_ALU_MODE_X_MUL_Y |
                   VECINT_BROADCAST_MODE(3) | VECINT_MULTIPLE_4);
#endif

#ifdef CPU_M1
        AMX_VECINT(vecint_args | VECINT_X_REG(0) | VECINT_Z_ROW(0) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
        AMX_VECINT(vecint_args | VECINT_X_REG(1) | VECINT_Z_ROW(16) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
        AMX_VECINT(vecint_args | VECINT_X_REG(2) | VECINT_Z_ROW(32) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
        AMX_VECINT(vecint_args | VECINT_X_REG(3) | VECINT_Z_ROW(48) | VECINT_ENABLE_MODE(0) | VECINT_ENABLE_VALUE(5));
#else
        AMX_VECINT(vecint_args | VECINT_X_REG(0) | VECINT_Z_ROW(0) | VECINT_BROADCAST_MODE(5) | VECINT_MULTIPLE_4);
#endif

#ifdef CPU_M1
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(0) | EXTRH_COPY_ONLY_Z_ROW(0));
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(1) | EXTRH_COPY_ONLY_Z_ROW(16));
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(2) | EXTRH_COPY_ONLY_Z_ROW(32));
        AMX_EXTRH(extrh_args | EXTRH_COPY_ONLY_REG(3) | EXTRH_COPY_ONLY_Z_ROW(48));
#else
        AMX_EXTRH(EXTRH_GENERIC | EXTRH_REG(0) | EXTRH_DESTINATION_X | EXTRH_Z_ROW(0) | EXTRH_LANE_WIDTH_16_BIT |
                  EXTRH_MULTIPLE_4);
#endif

        // Table lookup
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(0) | GENLUT_DST_REG(0) | GENLUT_TABLE_REG(0));
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(1) | GENLUT_DST_REG(1) | GENLUT_TABLE_REG(0));
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(2) | GENLUT_DST_REG(2) | GENLUT_TABLE_REG(0));
        AMX_GENLUT(genlut_generate_args | GENLUT_SRC_REG(3) | GENLUT_DST_REG(3) | GENLUT_TABLE_REG(0));

        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(0) | GENLUT_DST_REG(0) | GENLUT_TABLE_REG(1));
        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(1) | GENLUT_DST_REG(1) | GENLUT_TABLE_REG(1));
        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(2) | GENLUT_DST_REG(2) | GENLUT_TABLE_REG(1));
        AMX_GENLUT(genlut_lookup_args | GENLUT_SRC_REG(3) | GENLUT_DST_REG(3) | GENLUT_TABLE_REG(1));

        AMX_STX(AMX_PTR(&s[i]) | STX_REG(0) | STX_STORE_PAIR);
        AMX_STX(AMX_PTR(&s[i + 64]) | STX_REG(2) | STX_STORE_PAIR);
    }
}

void frodo_sample_n(uint16_t *s, const size_t n) {
    if (n & 127) {
        frodo_sample_n_64(s, n);
    }
    else {
        frodo_sample_n_128(s, n);
    }
}
