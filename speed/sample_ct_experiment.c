#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aarch64.h"
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

static int cmp_uint64(const void *a, const void *b) {
    return ((*((const uint64_t *)a)) - ((*((const uint64_t *)b))));
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

#if defined(AMX)
#include "noise_amx.c"
#elif defined(OPT_NEON)
#include "noise_neon.c"
#else
#include "noise.c"
#endif

// We use this function to initialize AMX. The use of __attribute__((constructor)) ensures this function is called
// before main().
__attribute__((constructor)) void init_AMX(void) {
    AMX_SET();
}

uint16_t *s_0, *s_random, *s_FFFE, *s_FFFF;

__attribute__((constructor)) void alloc_arrays(void) {
    s_0 = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
    s_random = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
    s_FFFE = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
    s_FFFF = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(uint16_t));
}

int main() {
    unsigned char entropy_input[48];

    SETUP_COUNTER();

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    randombytes((uint8_t *)s_random, PARAMS_N * 8 * sizeof(uint16_t));

    memset(s_0, 0, PARAMS_N * 8 * sizeof(uint16_t));

    for (int i = 0; i < PARAMS_N * 8; i++) {
        s_FFFE[i] = 0xFFFE;
    }

    memset(s_FFFF, 0xFFFF, PARAMS_N * 8 * sizeof(uint16_t));

    WRAP_FUNC("Cycles for random inputs: " CYCLE_TYPE " \n", cycles, time0, time1,
              frodo_sample_n(s_random, PARAMS_N * 8));

    WRAP_FUNC("Cycles for inputs = 0x0000: " CYCLE_TYPE " \n", cycles, time0, time1, frodo_sample_n(s_0, PARAMS_N * 8));

    WRAP_FUNC("Cycles for inputs = 0xFFFE: " CYCLE_TYPE " \n", cycles, time0, time1,
              frodo_sample_n(s_FFFE, PARAMS_N * 8));

    WRAP_FUNC("Cycles for inputs = 0xFFFF: " CYCLE_TYPE " \n", cycles, time0, time1,
              frodo_sample_n(s_FFFF, PARAMS_N * 8));

    return 0;
}
