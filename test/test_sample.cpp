#include <stddef.h>
#include <stdint.h>

#include "aarch64.h"
#include "gtest/gtest.h"
#include "test.h"

#if PARAMS_N == 640
// From PQCrypto-LWEKE/FrodoKEM/src/frodo640.c
static const uint16_t CDF_TABLE[13] = {4643,  13363, 20579, 25843, 29227, 31145, 32103,
                                       32525, 32689, 32745, 32762, 32766, 32767};
static const uint16_t CDF_TABLE_LEN = 13;
#elif PARAMS_N == 976
// From PQCrypto-LWEKE/FrodoKEM/src/frodo976.c
static const uint16_t CDF_TABLE[11] = {5638, 15915, 23689, 28571, 31116, 32217, 32613, 32731, 32760, 32766, 32767};
static const uint16_t CDF_TABLE_LEN = 11;
#elif PARAMS_N == 1344
// From PQCrypto-LWEKE/FrodoKEM/src/frodo1344.c
static const uint16_t CDF_TABLE[7] = {9142, 23462, 30338, 32361, 32725, 32765, 32767};
static const uint16_t CDF_TABLE_LEN = 7;
#else
#error Invalid/undefined PARAMS_N
#endif

#include "noise.c"

#define frodo_sample_n neon_frodo_sample_n
#include "noise_neon.c"
#undef frodo_sample_n

#define frodo_sample_n amx_frodo_sample_n
#include "noise_amx.c"
#undef frodo_sample_n

// We use this function to initialize AMX. The use of __attribute__((constructor)) ensures this function is called
// before main().
__attribute__((constructor)) void init_AMX(void) {
    AMX_SET();
}

static void test_all_cdf_table_values(void (*sample)(uint16_t *, const size_t)) {
    uint16_t s[65536] __attribute__((aligned(128))), s2[65536] __attribute__((aligned(128)));

    for (size_t i = 0; i < 65536; i++) {
        s2[i] = s[i] = i;
    }

    frodo_sample_n(s2, 65536);
    sample(s, 65536);

    EXPECT_TRUE(ArraysMatch(s2, s));
}

#define frodokem_sample_n3(IMPL, PARAMS_N) IMPL##_frodokem_##PARAMS_N##_sample_n
#define frodokem_sample_n2(IMPL, PARAMS_N) frodokem_sample_n3(IMPL, PARAMS_N)
#define frodokem_sample_n(IMPL) frodokem_sample_n2(IMPL, PARAMS_N)

TEST(frodokem_sample_n(amx), test_all_cdf_table_values) {
    test_all_cdf_table_values(amx_frodo_sample_n);
}

TEST(frodokem_sample_n(neon), test_all_cdf_table_values) {
    test_all_cdf_table_values(neon_frodo_sample_n);
}
