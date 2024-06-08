
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/aes/aes.h"
#include "config.h"
#include "fips202x2.h"
#include "frodo_macrify_4x.h"
#include "memory_alloc.h"
#include "rng.h"
#include "speed.h"

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

#undef __AVERAGE__
#define __MEDIAN__

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

#define PARAMS_STRIPE_STEP 8
#define PARAMS_NBAR 8
#define BYTES_SEED_A 16

#ifdef ALLOC_MMAP
int16_t *A_sa_4x;

__attribute__((constructor)) static void alloc_sa_4x(void) {
    A_sa_4x = MEMORY_ALLOC(PARAMS_N * PARAMS_N * sizeof(int16_t));
}
#endif

int frodo_mul_add_sa_plus_e_matmul_only_4x(uint16_t *out[4], const uint16_t *s[4],
                                           uint16_t *e[4]) {  // Generate-and-multiply: generate matrix A (N x N)
                                                              // column-wise, multiply by s' on the left. Inputs: s', e'
                                                              // (N_BAR x N) Output: out = s'*A + e' (N_BAR x N)
    int i, j, k, n;
#ifdef ALLOC_MMAP
    int16_t *A = A_sa_4x;
#else
    int16_t A[PARAMS_N * PARAMS_N];
#endif

    for (i = 0; i < PARAMS_N * PARAMS_N; i++) {
        A[i] = LE_TO_UINT16(A[i]);
    }

    for (n = 0; n < 4; n++) {
        memcpy(out[n], e[n], PARAMS_NBAR * PARAMS_N * sizeof(uint16_t));

        for (i = 0; i < PARAMS_N; i++) {  // Matrix multiplication-addition A*s + e
            for (k = 0; k < PARAMS_NBAR; k++) {
                uint16_t sum = 0;
                for (j = 0; j < PARAMS_N; j++) {
                    sum += A[j * PARAMS_N + i] * s[n][k * PARAMS_N + j];
                }
                out[n][k * PARAMS_N + i] += sum;  // Adding e. No need to reduce modulo 2^15, extra bits are taken care
                                                  // of during packing later on.
            }
        }
    }

    DoNotOptimize(out);

    return 1;
}

int main() {
    unsigned char entropy_input[48], seedA[BYTES_SEED_A];
#ifdef ALLOC_MMAP
    uint16_t *Bp[4], *Bpm[4], *Sp[4], *Spm[4], *Ep[4], *Epm[4];
    for (int i = 0; i < 4; i++) {
        Bp[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
        Bpm[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
        Sp[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
        Spm[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
        Ep[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
        Epm[i] = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    }
#else
    uint16_t Bp0[PARAMS_N * PARAMS_NBAR];
    uint16_t Bp1[PARAMS_N * PARAMS_NBAR];
    uint16_t Bp2[PARAMS_N * PARAMS_NBAR];
    uint16_t Bp3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Bp[4] = {Bp0, Bp1, Bp2, Bp3};

    uint16_t Bpm0[PARAMS_N * PARAMS_NBAR];
    uint16_t Bpm1[PARAMS_N * PARAMS_NBAR];
    uint16_t Bpm2[PARAMS_N * PARAMS_NBAR];
    uint16_t Bpm3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Bpm[4] = {Bpm0, Bpm1, Bpm2, Bpm3};

    uint16_t Sp0[PARAMS_N * PARAMS_NBAR];
    uint16_t Sp1[PARAMS_N * PARAMS_NBAR];
    uint16_t Sp2[PARAMS_N * PARAMS_NBAR];
    uint16_t Sp3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Sp[4] = {Sp0, Sp1, Sp2, Sp3};

    uint16_t Spm0[PARAMS_N * PARAMS_NBAR];
    uint16_t Spm1[PARAMS_N * PARAMS_NBAR];
    uint16_t Spm2[PARAMS_N * PARAMS_NBAR];
    uint16_t Spm3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Spm[4] = {Spm0, Spm1, Spm2, Spm3};

    uint16_t Ep0[PARAMS_N * PARAMS_NBAR];
    uint16_t Ep1[PARAMS_N * PARAMS_NBAR];
    uint16_t Ep2[PARAMS_N * PARAMS_NBAR];
    uint16_t Ep3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Ep[4] = {Ep0, Ep1, Ep2, Ep3};

    uint16_t Epm0[PARAMS_N * PARAMS_NBAR];
    uint16_t Epm1[PARAMS_N * PARAMS_NBAR];
    uint16_t Epm2[PARAMS_N * PARAMS_NBAR];
    uint16_t Epm3[PARAMS_N * PARAMS_NBAR];

    uint16_t *Epm[4] = {Epm0, Epm1, Epm2, Epm3};
#endif

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);
    randombytes(seedA, BYTES_SEED_A);

    SETUP_COUNTER();

    WRAP_FUNC("FrodoKEM s*A + e 4x: " CYCLE_TYPE "\n", cycles, time0, time1,
              frodo_mul_add_sa_plus_e_4x(Bp, (const uint16_t **)Sp, Ep, seedA));
    WRAP_FUNC("FrodoKEM s*A + e 4x (matmul only): " CYCLE_TYPE "\n", cycles, time0, time1,
              frodo_mul_add_sa_plus_e_matmul_only_4x(Bpm, (const uint16_t **)Spm, Epm));

    return 0;
}
