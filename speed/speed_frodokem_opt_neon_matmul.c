
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/aes/aes.h"
#include "config.h"
#include "fips202x2.h"
#include "frodo_macrify.h"
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

__attribute__((constructor)) static void alloc_as_speed(void) {
    a_row_as_speed = MEMORY_ALLOC(4 * PARAMS_N * sizeof(int16_t));
}
#endif

int frodo_mul_add_as_plus_e_matmul_only(uint16_t *out, const uint16_t *s,
                                        const uint16_t *e) {  // Generate-and-multiply: generate matrix A (N x N)
                                                              // row-wise, multiply by s on the right. Inputs: s, e (N x
                                                              // N_BAR) Output: out = A*s + e (N x N_BAR)
    int i, j, k;
#ifdef ALLOC_MMAP
    int16_t *a_row = a_row_as_speed;
#else
    ALIGN_HEADER(32) int16_t a_row[4 * PARAMS_N] ALIGN_FOOTER(32) = {0};
#endif

#ifndef OPT_NEON
    for (i = 0; i < (PARAMS_N * PARAMS_NBAR); i += 2) {
        *((uint32_t *)&out[i]) = *((uint32_t *)&e[i]);
    }
#endif

    for (i = 0; i < PARAMS_N; i += 4) {
#ifdef OPT_NEON
        uint16x8_t vA[2], vs[8], vout[2][8], ve[2];
        int p;

        for (k = 0; k < 4; k += 2) {
            for (j = 0; j < 8; j++) {
                vout[0][j] = vdupq_n_u16(0);
                vout[1][j] = vdupq_n_u16(0);
            }

            for (j = 0; j < PARAMS_N; j += 8) {
                vA[0] = vld1q_u16((uint16_t *)&a_row[k * PARAMS_N + j]);
                vA[1] = vld1q_u16((uint16_t *)&a_row[(k + 1) * PARAMS_N + j]);

                for (p = 0; p < 8; p++) {
                    vs[p] = vld1q_u16(&s[p * PARAMS_N + j]);
                    vout[0][p] = vmlaq_u16(vout[0][p], vA[0], vs[p]);
                    vout[1][p] = vmlaq_u16(vout[1][p], vA[1], vs[p]);
                }
            }

            for (j = 0; j < 2; j++) {
                vout[j][0] = vpaddq_u16(vout[j][0], vout[j][1]);
                vout[j][2] = vpaddq_u16(vout[j][2], vout[j][3]);
                vout[j][4] = vpaddq_u16(vout[j][4], vout[j][5]);
                vout[j][6] = vpaddq_u16(vout[j][6], vout[j][7]);

                vout[j][0] = vpaddq_u16(vout[j][0], vout[j][2]);
                vout[j][4] = vpaddq_u16(vout[j][4], vout[j][6]);

                vout[j][0] = vpaddq_u16(vout[j][0], vout[j][4]);

                ve[j] = vld1q_u16(&e[(i + k + j) * PARAMS_NBAR]);

                vout[j][0] = vaddq_u16(vout[j][0], ve[j]);

                vst1q_u16(&out[(i + k + j) * PARAMS_NBAR], vout[j][0]);
            }
        }
    }
#else
        for (k = 0; k < 4 * PARAMS_N; k++) {
            a_row[k] = LE_TO_UINT16(a_row[k]);
        }
        for (k = 0; k < PARAMS_NBAR; k++) {
            uint16_t sum[4] = {0};
            for (j = 0; j < PARAMS_N; j++) {  // Matrix-vector multiplication
                uint16_t sp = s[k * PARAMS_N + j];
                sum[0] += a_row[0 * PARAMS_N + j] * sp;  // Go through four lines with same s
                sum[1] += a_row[1 * PARAMS_N + j] * sp;
                sum[2] += a_row[2 * PARAMS_N + j] * sp;
                sum[3] += a_row[3 * PARAMS_N + j] * sp;
            }
            out[(i + 0) * PARAMS_NBAR + k] += sum[0];
            out[(i + 2) * PARAMS_NBAR + k] += sum[2];
            out[(i + 1) * PARAMS_NBAR + k] += sum[1];
            out[(i + 3) * PARAMS_NBAR + k] += sum[3];
        }
    }
#endif

    DoNotOptimize(out);

    return 1;
}

#ifdef ALLOC_MMAP
uint16_t *A_sa_speed;

__attribute__((constructor)) static void alloc_sa_speed(void) {
    A_sa_speed = MEMORY_ALLOC(PARAMS_N * 8 * sizeof(int16_t));
}
#endif

int frodo_mul_add_sa_plus_e_matmul_only(
    uint16_t *out, const uint16_t *s,
    uint16_t *e) {  // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
                    // Inputs: s', e' (N_BAR x N)
                    // Output: out = s'*A + e' (N_BAR x N)
                    // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in:
                    // J.W. Bos, M. Ofner, J. Renes, T. Schneider, C. van Vredendaal, "The Matrix Reloaded:
                    // Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, q, p;
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_speed;
#else
    ALIGN_HEADER(32) uint16_t A[PARAMS_N * 8] ALIGN_FOOTER(32);
#endif

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i += 8) {
#if defined(OPT_NEON)
        uint16x8_t vsp[8];

        for (j = 0; j < PARAMS_NBAR; j++) {
            vsp[j] = vld1q_u16(&s[j * PARAMS_N + i]);
        }

        for (q = 0; q < PARAMS_N; q += 8) {
            uint16x8_t vsum[8], vA[8];

            for (p = 0; p < 8; p++) {
                vA[p] = vld1q_u16(&A[p * PARAMS_N + q]);
            }

            for (j = 0; j < PARAMS_NBAR; j++) {
                vsum[j] = vld1q_u16(&e[j * PARAMS_N + q]);

                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[0], vsp[j], 0);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[1], vsp[j], 1);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[2], vsp[j], 2);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[3], vsp[j], 3);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[4], vsp[j], 4);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[5], vsp[j], 5);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[6], vsp[j], 6);
                vsum[j] = vmlaq_laneq_u16(vsum[j], vA[7], vsp[j], 7);

                vst1q_u16(&e[j * PARAMS_N + q], vsum[j]);
            }
        }
    }
#elif !defined(USE_AVX2)
        for (j = 0; j < PARAMS_NBAR; j++) {
            uint16_t sum = 0;
            int16_t sp[8];
            for (p = 0; p < 8; p++) {
                sp[p] = s[j * PARAMS_N + i + p];
            }
            for (q = 0; q < PARAMS_N; q++) {
                sum = e[j * PARAMS_N + q];
                for (p = 0; p < 8; p++) {
                    sum += sp[p] * A[p * PARAMS_N + q];
                }
                e[j * PARAMS_N + q] = sum;
            }
        }
    }
#endif
    memcpy((unsigned char *)out, (unsigned char *)e, 2 * PARAMS_N * PARAMS_NBAR);

    DoNotOptimize(out);

    return 1;
}

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
