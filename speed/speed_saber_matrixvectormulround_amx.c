#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SABER_params.h"
#include "aarch64.h"
#include "poly.h"
#include "polymodmul.h"
#include "rng.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

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

#ifdef ALLOC_MMAP
static uint16_t *A_, *s_, *b_;

__attribute__((constructor)) static void alloc_kp(void) {
    A_ = MEMORY_ALLOC(SABER_L * SABER_L * SABER_N * sizeof(uint16_t));

    s_ = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    b_ = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
}
#endif

int main() {
#ifdef ALLOC_STACK
    uint16_t A[SABER_L][SABER_L][SABER_N];

    uint16_t s[SABER_L][SABER_N];
    uint16_t b[SABER_L][SABER_N] = {0};
#else
    uint16_t(*A)[SABER_L][SABER_N] = (uint16_t(*)[SABER_L][SABER_N])A_;

    uint16_t(*s)[SABER_N] = (uint16_t(*)[SABER_N])s_;
    uint16_t(*b)[SABER_N] = (uint16_t(*)[SABER_N])b_;

    memset(b, 0, SABER_L * SABER_N * sizeof(uint16_t));
#endif

    unsigned char entropy_input[48] = {0};

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);
    // randombytes((unsigned char *)x, 256 * sizeof(uint16_t));
    // randombytes((unsigned char *)y, 256 * sizeof(uint16_t));

    SETUP_COUNTER();

    WRAP_FUNC("MatrixVectorMulRound: " CYCLE_TYPE "\n", cycles, time0, time1,
              MatrixVectorMulRound((const uint16_t(*)[SABER_L][SABER_N])A, (const uint16_t(*)[SABER_N])s, b, 1));

    return 0;
}
