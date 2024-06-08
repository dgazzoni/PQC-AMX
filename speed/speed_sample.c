
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "memory_alloc.h"
#include "rng.h"

#ifndef NTESTS
#error Please #define NTESTS to the desired number of repetitions for the benchmark
#endif

uint64_t time0, time1;
uint64_t cycles[NTESTS];

#ifdef __APPLE__

#include "m1cycles.h"
#define SETUP_COUNTER() \
    {                   \
        (void)cycles;   \
        setup_rdtsc();  \
    }
#define CYCLE_TYPE "%lld"
#define GET_TIME rdtsc()

#else

#include "hal.h"
#define SETUP_COUNTER() \
    {}
#define CYCLE_TYPE "%ld"
#define GET_TIME hal_get_time()

#endif

#undef __MEDIAN__
#define __AVERAGE__  // For these smaller tests, especially those involving AMX, medians have proved far too optimistic

#ifdef __AVERAGE__

#define LOOP_INIT(__clock0, __clock1) \
    { __clock0 = GET_TIME; }
#define LOOP_TAIL(__f_string, records, __clock0, __clock1)  \
    {                                                       \
        __clock1 = GET_TIME;                                \
        printf(__f_string, (__clock1 - __clock0) / NTESTS); \
    }
#define BODY_INIT(__clock0, __clock1) \
    {}
#define BODY_TAIL(records, __clock0, __clock1) \
    {}

#elif defined(__MEDIAN__)

static int cmp_uint64(const void *a, const void *b) {
    return ((*((const uint64_t *)a)) - ((*((const uint64_t *)b))));
}

#define LOOP_INIT(__clock0, __clock1) \
    {}
#define LOOP_TAIL(__f_string, records, __clock0, __clock1)    \
    {                                                         \
        qsort(records, sizeof(uint64_t), NTESTS, cmp_uint64); \
        printf(__f_string, records[NTESTS >> 1]);             \
    }
#define BODY_INIT(__clock0, __clock1) \
    { __clock0 = GET_TIME; }
#define BODY_TAIL(records, __clock0, __clock1) \
    {                                          \
        __clock1 = GET_TIME;                   \
        records[i] = __clock1 - __clock0;      \
    }

#endif

#define WRAP_FUNC(__f_string, records, __clock0, __clock1, func) \
    {                                                            \
        LOOP_INIT(__clock0, __clock1);                           \
        for (size_t i = 0; i < NTESTS; i++) {                    \
            BODY_INIT(__clock0, __clock1);                       \
            func;                                                \
            BODY_TAIL(records, __clock0, __clock1);              \
        }                                                        \
        LOOP_TAIL(__f_string, records, __clock0, __clock1);      \
    }

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

int main() {
    unsigned char entropy_input[48];
    uint16_t *s_opt = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
    uint16_t *s_NEON = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
    uint16_t *s_AMX = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);
    randombytes((unsigned char *)s_opt, PARAMS_N * 8 * sizeof(uint16_t));
    randombytes((unsigned char *)s_NEON, PARAMS_N * 8 * sizeof(uint16_t));
    randombytes((unsigned char *)s_AMX, PARAMS_N * 8 * sizeof(uint16_t));

    SETUP_COUNTER();

    WRAP_FUNC("FrodoKEM sample: " CYCLE_TYPE "\n", cycles, time0, time1, frodo_sample_n(s_opt, PARAMS_N * 8));
    WRAP_FUNC("FrodoKEM neon_sample: " CYCLE_TYPE "\n", cycles, time0, time1,
              neon_frodo_sample_n(s_NEON, PARAMS_N * 8));
    WRAP_FUNC("FrodoKEM amx_sample: " CYCLE_TYPE "\n", cycles, time0, time1, amx_frodo_sample_n(s_AMX, PARAMS_N * 8));

    return 0;
}
