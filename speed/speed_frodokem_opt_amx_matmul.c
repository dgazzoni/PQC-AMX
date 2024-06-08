
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/aes/aes.h"
#include "config.h"
#include "fips202x2.h"
#include "frodo_macrify.h"
#include "frodokem_opt_matmul.h"
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
int16_t *a_row_as_speed;
uint16_t *s_as_speed;

__attribute__((constructor)) static void alloc_as_speed(void) {
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    a_row_as_speed = MEMORY_ALLOC(32 * n32 * sizeof(int16_t));
    s_as_speed = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_as_plus_e_matmul_only(uint16_t *out, const uint16_t *sT,
                                        const uint16_t *e) {  // Generate-and-multiply: generate matrix A (N x N)
                                                              // row-wise, multiply by s on the right. Inputs: s, e (N x
                                                              // N_BAR) Output: out = A*s + e (N x N_BAR)
    int i, j;

#ifdef ALLOC_MMAP
    int16_t *a_row = a_row_as_speed;
    uint16_t *s = s_as_speed;
#else
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    ALIGN_HEADER(32) int16_t a_row[32 * n32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif

    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
        }
    }

    for (i = 0; i < PARAMS_N; i += 32) {
        amx_mul_add_as_plus_e_up_to_32_rows(out, a_row, s, e, PARAMS_N, i);
    }

    return 1;
}

#if defined(USE_AES128_FOR_A)
#ifdef ALLOC_MMAP
uint16_t *A_sa_speed, *s_sa_speed;

__attribute__((constructor)) static void alloc_sa_speed(void) {
    A_sa_speed = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(int16_t));
    s_sa_speed = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e_matmul_only(
    uint16_t *out, const uint16_t *sT,
    uint16_t *e) {  // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
                    // Inputs: s', e' (N_BAR x N)
                    // Output: out = s'*A + e' (N_BAR x N)
                    // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in:
                    // J.W. Bos, M. Ofner, J. Renes, T. Schneider, C. van Vredendaal, "The Matrix Reloaded:
                    // Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j;
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_speed;
    uint16_t *s = s_sa_speed;
#else
    ALIGN_HEADER(32) uint16_t A[PARAMS_N * 32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif

    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
        }
    }

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i += 32) {
        amx_mul_add_sa_plus_e_up_to_32_cols(out, (int16_t *)A, s, e, PARAMS_N, i);
    }

    DoNotOptimize(out);

    return 1;
}
#elif defined(USE_SHAKE128_FOR_A)  // SHAKE128
#ifdef ALLOC_MMAP
uint16_t *A_sa_speed, *s_sa_speed;

__attribute__((constructor)) static void alloc_sa_speed(void) {
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    A_sa_speed = MEMORY_ALLOC(32 * n32 * sizeof(int16_t));
    s_sa_speed = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e_matmul_only(
    uint16_t *out, const uint16_t *sT,
    uint16_t *e) {  // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
                    // Inputs: s', e' (N_BAR x N)
                    // Output: out = s'*A + e' (N_BAR x N)
                    // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in:
                    // J.W. Bos, M. Ofner, J. Renes, T. Schneider, C. van Vredendaal, "The Matrix Reloaded:
                    // Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j;
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_speed;
    uint16_t *s = s_sa_speed;
#else
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    ALIGN_HEADER(32) uint16_t A[32 * n32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif

    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
        }
    }

    for (i = 0; i < PARAMS_N; i += 32) {
        amx_mul_add_sa_plus_e_up_to_32_rows(out, (int16_t *)A, s, e, PARAMS_N, i);
    }

    DoNotOptimize(out);

    return 1;
}
#endif

int main() {
    unsigned char entropy_input[48], seedA[BYTES_SEED_A];
#ifdef ALLOC_MMAP
    uint16_t *B = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Bm = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Bp = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Bpm = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *S = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Sm = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Sp = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Spm = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *E = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Em = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Ep = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
    uint16_t *Epm = MEMORY_ALLOC(PARAMS_N * PARAMS_NBAR * sizeof(uint16_t));
#else
    uint16_t B[PARAMS_N * PARAMS_NBAR];
    uint16_t Bm[PARAMS_N * PARAMS_NBAR];
    uint16_t Bp[PARAMS_N * PARAMS_NBAR];
    uint16_t Bpm[PARAMS_N * PARAMS_NBAR];
    uint16_t S[PARAMS_N * PARAMS_NBAR];
    uint16_t Sm[PARAMS_N * PARAMS_NBAR];
    uint16_t Sp[PARAMS_N * PARAMS_NBAR];
    uint16_t Spm[PARAMS_N * PARAMS_NBAR];
    uint16_t E[PARAMS_N * PARAMS_NBAR];
    uint16_t Em[PARAMS_N * PARAMS_NBAR];
    uint16_t Ep[PARAMS_N * PARAMS_NBAR];
    uint16_t Epm[PARAMS_N * PARAMS_NBAR];
#endif

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);
    randombytes(seedA, BYTES_SEED_A);

    SETUP_COUNTER();

    WRAP_FUNC("FrodoKEM A*s + e: " CYCLE_TYPE "\n", cycles, time0, time1, frodo_mul_add_as_plus_e(B, S, E, seedA));
    WRAP_FUNC("FrodoKEM A*s + e (matmul only): " CYCLE_TYPE "\n", cycles, time0, time1,
              frodo_mul_add_as_plus_e_matmul_only(Bm, Sm, Em));
    WRAP_FUNC("FrodoKEM s*A + e: " CYCLE_TYPE "\n", cycles, time0, time1, frodo_mul_add_sa_plus_e(Bp, Sp, Ep, seedA));
    WRAP_FUNC("FrodoKEM s*A + e (matmul only): " CYCLE_TYPE "\n", cycles, time0, time1,
              frodo_mul_add_sa_plus_e_matmul_only(Bpm, Spm, Epm));

    return 0;
}
