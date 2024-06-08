/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: matrix arithmetic functions used by the KEM
*********************************************************************************************/

#if defined(USE_AES128_FOR_A)
#if !defined(USE_OPENSSL)
    #include "../../common/aes/aes.h"
#else
    #include "../../common/aes/aes_openssl.h"
#endif
#elif defined (USE_SHAKE128_FOR_A)
#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || (__ARM_FEATURE_SHA3)
    #include "fips202x2.h"
#elif !defined(USE_AVX2)
    #include "../../common/sha3/fips202.h"
#else
    #include "../../common/sha3/fips202x4.h"
#endif
#endif    
#if defined(USE_AVX2)
    #include <immintrin.h>
#endif

#ifdef OPT_NEON
void AES128_ECB_genA_rows(const uint8_t *schedule, uint16_t *A, int n_rows, int first_row, int n, int n_cols);
#endif

#ifdef ALLOC_MMAP
uint16_t *A_sa_4x;

__attribute__((constructor)) static void alloc_sa_4x(void) {
    A_sa_4x = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(int16_t));
}
#endif

int frodo_mul_add_sa_plus_e_4x(uint16_t *out[4], const uint16_t *s[4], uint16_t *e[4], const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, q, p, n; 
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa_4x;
#else
    ALIGN_HEADER(32) uint16_t A[PARAMS_N*8] ALIGN_FOOTER(32);
#endif
    memset(A, 0, PARAMS_N * 8 * sizeof(uint16_t));

#if defined(USE_AES128_FOR_A)
#if !defined(USE_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);
#else
    EVP_CIPHER_CTX *aes_key_schedule;
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();
#endif
    // Initialize matrix used for encryption
#ifndef OPT_NEON
    ALIGN_HEADER(32) uint16_t Ainit[PARAMS_N*8] ALIGN_FOOTER(32) = {0};
       
    for(j = 0; j < PARAMS_N; j+=8) {
        Ainit[0*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[1*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[2*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[3*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[4*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[5*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[6*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        Ainit[7*PARAMS_N + j + 1] = UINT16_TO_LE(j);
    }
#endif

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
#ifdef OPT_NEON
        AES128_ECB_genA_rows(aes_key_schedule, A, 8, i, PARAMS_N, PARAMS_N);
#else
        // Generate 8 rows of A on-the-fly using AES
        for (q = 0; q < 8; q++) {
            for (p = 0; p < PARAMS_N; p+=8) {
                Ainit[q*PARAMS_N + p] = UINT16_TO_LE(i+q);
            }
        }
        size_t A_len = 8 * PARAMS_N * sizeof(uint16_t);
#if !defined(USE_OPENSSL)
        AES128_ECB_enc_sch((uint8_t*)Ainit, A_len, aes_key_schedule, (uint8_t*)A);
#else   
        if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)A, &len, (uint8_t*)Ainit, A_len)) handleErrors();
#endif 
#endif
#elif defined (USE_SHAKE128_FOR_A)  // SHAKE128
#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || (__ARM_FEATURE_SHA3)
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t* seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 8) {
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        shake128x2((uint8_t*)(A + 0*PARAMS_N), (uint8_t*)(A + 1*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        seed_A_origin_0[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 3);
        shake128x2((uint8_t*)(A + 2*PARAMS_N), (uint8_t*)(A + 3*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        seed_A_origin_0[0] = UINT16_TO_LE(i + 4);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 5);
        shake128x2((uint8_t*)(A + 4*PARAMS_N), (uint8_t*)(A + 5*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        seed_A_origin_0[0] = UINT16_TO_LE(i + 6);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 7);
        shake128x2((uint8_t*)(A + 6*PARAMS_N), (uint8_t*)(A + 7*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);

#elif !defined(USE_AVX2)
    uint8_t seed_A_separated[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
    memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
        seed_A_origin[0] = UINT16_TO_LE(i + 0);
        shake128((unsigned char*)(A + 0*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 1);
        shake128((unsigned char*)(A + 1*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 2);
        shake128((unsigned char*)(A + 2*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 3);
        shake128((unsigned char*)(A + 3*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 4);
        shake128((unsigned char*)(A + 4*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 5);
        shake128((unsigned char*)(A + 5*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 6);
        shake128((unsigned char*)(A + 6*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 7);
        shake128((unsigned char*)(A + 7*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A); 
#else  // Using vector intrinsics
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_2[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_3[2 + BYTES_SEED_A];
    uint16_t *seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t *seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    uint16_t *seed_A_origin_2 = (uint16_t*)&seed_A_separated_2;
    uint16_t *seed_A_origin_3 = (uint16_t*)&seed_A_separated_3;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_2[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_3[2], seed_A, BYTES_SEED_A);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i+=8) {
        // Generate hash output
        // First 4 rows
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 3);
        shake128_4x((unsigned char*)(A + 0*PARAMS_N), (unsigned char*)(A + 1*PARAMS_N), (unsigned char*)(A + 2*PARAMS_N), (unsigned char*)(A + 3*PARAMS_N),
                    (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
        // Second 4 rows
        seed_A_origin_0[0] = UINT16_TO_LE(i + 4);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 5);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 6);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 7);
        shake128_4x((unsigned char*)(A + 4*PARAMS_N), (unsigned char*)(A + 5*PARAMS_N), (unsigned char*)(A + 6*PARAMS_N), (unsigned char*)(A + 7*PARAMS_N),
                    (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
#endif
#endif

#if defined(OPT_NEON)
        for (n = 0; n < 4; n++) {
            uint16x8_t vsp[8];

            for (j = 0; j < PARAMS_NBAR; j++) {
                vsp[j] = vld1q_u16(&s[n][j * PARAMS_N + i]);
            }

            for (q = 0; q < PARAMS_N; q += 8) {
                uint16x8_t vsum[8], vA[8];

                for (p = 0; p < 8; p++) {
                    vA[p] = vld1q_u16(&A[p * PARAMS_N + q]);
                }

                for (j = 0; j < PARAMS_NBAR; j++) {
                    vsum[j] = vld1q_u16(&e[n][j * PARAMS_N + q]);

                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[0], vsp[j], 0);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[1], vsp[j], 1);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[2], vsp[j], 2);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[3], vsp[j], 3);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[4], vsp[j], 4);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[5], vsp[j], 5);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[6], vsp[j], 6);
                    vsum[j] = vmlaq_laneq_u16(vsum[j], vA[7], vsp[j], 7);

                    vst1q_u16(&e[n][j * PARAMS_N + q], vsum[j]);
                }
            }
        }
    }
#elif !defined(USE_AVX2)
        for (j = 0; j < PARAMS_NBAR; j++) {
            for (n = 0; n < 4; n++) {
                uint16_t sum = 0;
                int16_t sp[8];
                for (p = 0; p < 8; p++) {
                    sp[p] = s[n][j*PARAMS_N + i + p];
                }
                for (q = 0; q < PARAMS_N; q++) {
                    sum = e[n][j*PARAMS_N + q];
                    for (p = 0; p < 8; p++) {
                        sum += sp[p] * A[p*PARAMS_N + q];
                    }
                    e[n][j*PARAMS_N + q] = sum;
                }
            }
        }
    }
#endif

    for (n = 0; n < 4; n++) {
        memcpy((unsigned char*)out[n], (unsigned char*)e[n], 2*PARAMS_N*PARAMS_NBAR);
    }

#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}

void frodo_mul_add_sb_plus_e_4x(uint16_t *out[4], const uint16_t *b, const uint16_t *s[4], const uint16_t *e[4]) 
{
    for (int n = 0; n < 4; n++) {
        frodo_mul_add_sb_plus_e(out[n], b, s[n], e[n]);
    }
}
