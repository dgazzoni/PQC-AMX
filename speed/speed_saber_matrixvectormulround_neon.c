#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "NTT.h"
#include "SABER_params.h"
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
static uint32_t *A_NTT_, *s_NTT_asymmetric_, *s_NTT_;
static uint16_t *b_;

__attribute__((constructor)) static void alloc_kp(void) {
    A_NTT_ = MEMORY_ALLOC(SABER_L * SABER_L * SABER_N * sizeof(uint32_t));
    s_NTT_asymmetric_ = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint32_t));
    s_NTT_ = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    b_ = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
}
#endif

extern void __asm_round(uint16_t des[SABER_N], uint32_t src[SABER_N]);

void MatrixVectorMulRound(uint32_t A_NTT[SABER_L][SABER_L][SABER_N], uint32_t s_NTT_asymmetric[SABER_L][SABER_N],
                          uint32_t s_NTT[SABER_L][SABER_N], uint16_t b[SABER_L][SABER_N]) {
    for (int j = 0; j < SABER_L; j++) {
        for (int k = 0; k < SABER_L; k++) {
            NTT(A_NTT[j][k]);
        }
    }
    for (int j = 0; j < SABER_L; j++) {
        NTT_heavy(s_NTT_asymmetric[j], s_NTT[j]);
    }
    for (int j = 0; j < SABER_L; j++) {
        __asm_asymmetric_mul(&(A_NTT[j][0][0]), &(s_NTT[0][0]), &(s_NTT_asymmetric[0][0]), constants);
    }
    for (int j = 0; j < SABER_L; j++) {
        iNTT(&(A_NTT[j][0][0]));
    }
    for (int i = 0; i < SABER_L; i++) {
        __asm_round(b[i], A_NTT[i][0]);
    }
}

int main() {
#ifdef ALLOC_STACK
    uint32_t A_NTT[SABER_L][SABER_L][SABER_N];
    uint32_t s_NTT_asymmetric[SABER_L][SABER_N];
    uint32_t s_NTT[SABER_L][SABER_N];
    uint16_t b[SABER_L][SABER_N] = {0};
#else
    uint32_t(*A_NTT)[SABER_L][SABER_N] = (uint32_t(*)[SABER_L][SABER_N])A_NTT_;
    uint32_t(*s_NTT_asymmetric)[SABER_N] = (uint32_t(*)[SABER_N])s_NTT_asymmetric_;
    uint32_t(*s_NTT)[SABER_N] = (uint32_t(*)[SABER_N])s_NTT_;
    uint16_t(*b)[SABER_N] = (uint16_t(*)[SABER_N])b_;
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
              MatrixVectorMulRound(A_NTT, s_NTT_asymmetric, s_NTT, b));

    return 0;
}
