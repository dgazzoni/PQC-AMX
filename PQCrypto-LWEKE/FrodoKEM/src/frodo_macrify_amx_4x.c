/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: matrix arithmetic functions used by the KEM
*********************************************************************************************/

#if defined(USE_AES128_FOR_A)
    #include "../../common/aes/aes.h"
#elif defined (USE_SHAKE128_FOR_A)
    #include "fips202x2.h"
#endif
#include "frodokem_matmul_4x.h"
#include "frodokem_opt_matmul_4x.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

void AES128_ECB_genA_32cols(const uint8_t *schedule, uint16_t *A, int n, int first_col);

#if defined(USE_AES128_FOR_A)
#ifdef ALLOC_MMAP
uint16_t *A_sa_4x, *s_sa_4x;

__attribute__((constructor)) static void alloc_sa_4x(void) {
    A_sa_4x = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(int16_t));
    s_sa_4x = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e_4x(uint16_t *out[4], const uint16_t *sT[4], uint16_t *e[4], const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, n;
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_4x;
    uint16_t *s = s_sa_4x;
#else
    ALIGN_HEADER(32) uint16_t A[PARAMS_N * 32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[PARAMS_N * 32] ALIGN_FOOTER(32);
#endif
    memset(A, 0, PARAMS_N * 32 * sizeof(uint16_t));

    for (n = 0; n < 4; n++) {
        for (i = 0; i < PARAMS_N; i++) {
            for (j = 0; j < 8; j++) {
                s[i * 32 + (8*n + j)] = sT[n][j * PARAMS_N + i];
            }
        }
    }

    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i += 32) {
        AES128_ECB_genA_32cols(aes_key_schedule, A, i, PARAMS_N);

        amx_mul_add_sa_plus_e_up_to_32_cols_4x(out, (int16_t *)A, s, (const uint16_t**)e, PARAMS_N, i);
    }

    AES128_free_schedule(aes_key_schedule);

    return 1;
}
#elif defined (USE_SHAKE128_FOR_A)  // SHAKE128
#ifdef ALLOC_MMAP
uint16_t *A_sa_4x, *s_sa_4x;

__attribute__((constructor)) static void alloc_sa_4x(void) {
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    A_sa_4x = MEMORY_ALLOC(32 * n32 * sizeof(int16_t));
    s_sa_4x = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e_4x(uint16_t *out[4], const uint16_t *sT[4], uint16_t *e[4], const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, n;
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_4x;
    uint16_t *s = s_sa_4x;
#else
    ALIGN_HEADER(32) uint16_t A[32 * n32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[PARAMS_N * 32] ALIGN_FOOTER(32);
#endif
    memset(A, 0, 32 * n32 * sizeof(uint16_t));

    for (n = 0; n < 4; n++) {
        for (i = 0; i < PARAMS_N; i++) {
            for (j = 0; j < 8; j++) {
                s[i * 32 + (8*n + j)] = sT[n][j * PARAMS_N + i];
            }
        }
    }

    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t* seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 32) {
        for (j = 0; j < MIN(PARAMS_N - i, 32); j += 2) {
            seed_A_origin_0[0] = i + j;
            seed_A_origin_1[0] = i + j + 1;
            shake128x2((uint8_t*)(A + j*n32), (uint8_t*)(A + (j + 1)*n32), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        }

        amx_mul_add_sa_plus_e_up_to_32_rows_4x(out, (int16_t *)A, s, (const uint16_t**)e, PARAMS_N, i);
    }

    return 1;
}
#endif

void frodo_mul_add_sb_plus_e_4x(uint16_t *out[4], const uint16_t *b, const uint16_t *s[4], const uint16_t *e[4]) 
{
    amx_mul_add_sb_plus_e_4x(out, b, s, e, PARAMS_N);
}
