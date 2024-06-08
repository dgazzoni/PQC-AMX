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
int16_t *a_row_as, *a_row_temp_as;

__attribute__((constructor)) static void alloc_as(void) {
    a_row_as = MEMORY_ALLOC(4 * PARAMS_N * sizeof(int16_t));
    a_row_temp_as = MEMORY_ALLOC(4 * PARAMS_N * sizeof(int16_t));
}
#endif

int frodo_mul_add_as_plus_e(uint16_t *out, const uint16_t *s, const uint16_t *e, const uint8_t *seed_A) 
{ // Generate-and-multiply: generate matrix A (N x N) row-wise, multiply by s on the right.
  // Inputs: s, e (N x N_BAR)
  // Output: out = A*s + e (N x N_BAR)
    int i, j, k;
#ifdef ALLOC_MMAP
    int16_t *a_row = a_row_as;
#else
    ALIGN_HEADER(32) int16_t a_row[4*PARAMS_N] ALIGN_FOOTER(32) = {0};
#endif
    memset(a_row, 0, 4 * PARAMS_N * sizeof(int16_t));

#ifndef OPT_NEON
    for (i = 0; i < (PARAMS_N*PARAMS_NBAR); i += 2) {    
        *((uint32_t*)&out[i]) = *((uint32_t*)&e[i]);
    }    
#endif
    
#if defined(USE_AES128_FOR_A)
#ifndef OPT_NEON
#ifdef ALLOC_MMAP
    int16_t *a_row_temp = a_row_temp_as;
#else
    int16_t a_row_temp[4*PARAMS_N];                             // Take four lines of A at once 
#endif
    memset(a_row_temp, 0, 4 * PARAMS_N * sizeof(int16_t));
#endif
#if !defined(USE_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);   
#else
    EVP_CIPHER_CTX *aes_key_schedule;    
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();    
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();    
#endif

#ifndef OPT_NEON                                     
    for (j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {
        a_row_temp[j + 1 + 0*PARAMS_N] = UINT16_TO_LE(j);       // Loading values in the little-endian order
        a_row_temp[j + 1 + 1*PARAMS_N] = UINT16_TO_LE(j);
        a_row_temp[j + 1 + 2*PARAMS_N] = UINT16_TO_LE(j);
        a_row_temp[j + 1 + 3*PARAMS_N] = UINT16_TO_LE(j);
    }
#endif

    for (i = 0; i < PARAMS_N; i += 4) {
#ifdef OPT_NEON
        AES128_ECB_genA_rows(aes_key_schedule, (uint16_t*)a_row, 4, i, PARAMS_N, PARAMS_N);
#else
        for (j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {    // Go through A, four rows at a time
            a_row_temp[j + 0*PARAMS_N] = UINT16_TO_LE(i+0);     // Loading values in the little-endian order                                
            a_row_temp[j + 1*PARAMS_N] = UINT16_TO_LE(i+1);
            a_row_temp[j + 2*PARAMS_N] = UINT16_TO_LE(i+2);
            a_row_temp[j + 3*PARAMS_N] = UINT16_TO_LE(i+3);
        }

#if !defined(USE_OPENSSL)
        AES128_ECB_enc_sch((uint8_t*)a_row_temp, 4*PARAMS_N*sizeof(int16_t), aes_key_schedule, (uint8_t*)a_row);
#else   
        if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)a_row, &len, (uint8_t*)a_row_temp, 4*PARAMS_N*sizeof(int16_t))) handleErrors();
#endif
#endif
#elif defined (USE_SHAKE128_FOR_A)       
#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || (__ARM_FEATURE_SHA3)
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t* seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 4) {
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        shake128x2((uint8_t*)(a_row + 0*PARAMS_N), (uint8_t*)(a_row + 1*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        seed_A_origin_0[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 3);
        shake128x2((uint8_t*)(a_row + 2*PARAMS_N), (uint8_t*)(a_row + 3*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
#elif !defined(USE_AVX2)
    uint8_t seed_A_separated[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
    memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 4) {
        seed_A_origin[0] = UINT16_TO_LE(i + 0);
        shake128((unsigned char*)(a_row + 0*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 1);
        shake128((unsigned char*)(a_row + 1*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 2);
        shake128((unsigned char*)(a_row + 2*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
        seed_A_origin[0] = UINT16_TO_LE(i + 3);
        shake128((unsigned char*)(a_row + 3*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
#else
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_2[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_3[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t* seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    uint16_t* seed_A_origin_2 = (uint16_t*)&seed_A_separated_2;
    uint16_t* seed_A_origin_3 = (uint16_t*)&seed_A_separated_3;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_2[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_3[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i += 4) {
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        seed_A_origin_2[0] = UINT16_TO_LE(i + 2);
        seed_A_origin_3[0] = UINT16_TO_LE(i + 3);
        shake128_4x((unsigned char*)(a_row), (unsigned char*)(a_row + PARAMS_N), (unsigned char*)(a_row + 2*PARAMS_N), (unsigned char*)(a_row + 3*PARAMS_N), 
                    (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, seed_A_separated_2, seed_A_separated_3, 2 + BYTES_SEED_A);
#endif
#endif
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
            for (j = 0; j < PARAMS_N; j++) {                    // Matrix-vector multiplication            
                uint16_t sp = s[k*PARAMS_N + j];
                sum[0] += a_row[0*PARAMS_N + j] * sp;           // Go through four lines with same s
                sum[1] += a_row[1*PARAMS_N + j] * sp;
                sum[2] += a_row[2*PARAMS_N + j] * sp;
                sum[3] += a_row[3*PARAMS_N + j] * sp;
            }
            out[(i+0)*PARAMS_NBAR + k] += sum[0];
            out[(i+2)*PARAMS_NBAR + k] += sum[2];
            out[(i+1)*PARAMS_NBAR + k] += sum[1];
            out[(i+3)*PARAMS_NBAR + k] += sum[3];
        }
    }
#endif
    
#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}

#ifdef ALLOC_MMAP
uint16_t *A_sa;

__attribute__((constructor)) static void alloc_sa(void) {
    A_sa = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(int16_t));
}
#endif

int frodo_mul_add_sa_plus_e(uint16_t *out, const uint16_t *s, uint16_t *e, const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j, q, p; 
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa;
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
                sp[p] = s[j*PARAMS_N + i + p];
            }
            for (q = 0; q < PARAMS_N; q++) {
                sum = e[j*PARAMS_N + q];
                for (p = 0; p < 8; p++) {
                    sum += sp[p] * A[p*PARAMS_N + q];
                }
                e[j*PARAMS_N + q] = sum;
            }
        }
    }
#else  // Using vector intrinsics
        for (j = 0; j < PARAMS_NBAR; j++) {
            __m256i b, sp[8], acc;
            for (p = 0; p < 8; p++) {
                sp[p] = _mm256_set1_epi16(s[j*PARAMS_N + i + p]);
            }
            for (q = 0; q < PARAMS_N; q+=16) {
                acc = _mm256_load_si256((__m256i*)&e[j*PARAMS_N + q]);
                for (p = 0; p < 8; p++) {
                    b = _mm256_load_si256((__m256i*)&A[p*PARAMS_N + q]);
                    b = _mm256_mullo_epi16(b, sp[p]);
                    acc = _mm256_add_epi16(b, acc);
                }
                _mm256_store_si256((__m256i*)&e[j*PARAMS_N + q], acc);
            }
        }
    }
#endif
    memcpy((unsigned char*)out, (unsigned char*)e, 2*PARAMS_N*PARAMS_NBAR);

#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}


void frodo_mul_bs(uint16_t *out, const uint16_t *b, const uint16_t *s) 
{ // Multiply by s on the right
  // Inputs: b (N_BAR x N), s (N x N_BAR)
  // Output: out = b*s (N_BAR x N_BAR)
    int i, j, k;

#ifdef OPT_NEON
    int p;
    uint16x8_t vmask = vdupq_n_u16(PARAMS_Q - 1);

    for (i = 0; i < PARAMS_NBAR; i += 4) {
        uint16x8_t vb[2], vs[8], vout[2][8];

        for (k = 0; k < 4; k += 2) {
            for (j = 0; j < 8; j++) {
                vout[0][j] = vdupq_n_u16(0);
                vout[1][j] = vdupq_n_u16(0);
            }

            for (j = 0; j < PARAMS_N; j += 8) {
                vb[0] = vld1q_u16((uint16_t *)&b[(i + k) * PARAMS_N + j]);
                vb[1] = vld1q_u16((uint16_t *)&b[(i + k + 1) * PARAMS_N + j]);

                for (p = 0; p < 8; p++) {
                    vs[p] = vld1q_u16(&s[p * PARAMS_N + j]);
                    vout[0][p] = vmlaq_u16(vout[0][p], vb[0], vs[p]);
                    vout[1][p] = vmlaq_u16(vout[1][p], vb[1], vs[p]);
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

                vout[j][0] = vandq_u16(vout[j][0], vmask);

                vst1q_u16(&out[(i + k + j) * PARAMS_NBAR], vout[j][0]);
            }
        }
    }
#else
    for (i = 0; i < PARAMS_NBAR; i++) {
        for (j = 0; j < PARAMS_NBAR; j++) {
            out[i*PARAMS_NBAR + j] = 0;
            for (k = 0; k < PARAMS_N; k++) {
                out[i*PARAMS_NBAR + j] += b[i*PARAMS_N + k] * (int16_t)s[j*PARAMS_N + k];
            }
            out[i*PARAMS_NBAR + j] = (uint32_t)(out[i*PARAMS_NBAR + j]) & ((1<<PARAMS_LOGQ)-1);
        }
    }
#endif
}


void frodo_mul_add_sb_plus_e(uint16_t *out, const uint16_t *b, const uint16_t *s, const uint16_t *e) 
{ // Multiply by s on the left
  // Inputs: b (N x N_BAR), s (N_BAR x N), e (N_BAR x N_BAR)
  // Output: out = s*b + e (N_BAR x N_BAR)
    int i, j;

#ifdef OPT_NEON
    uint16x8_t vout[8], vs, vb[8], vmask = vdupq_n_u16(PARAMS_Q - 1);

    for (i = 0; i < 8; i++) {
        vout[i] = vld1q_u16(&e[i * PARAMS_NBAR]);
    }

    for (j = 0; j < PARAMS_N; j += 8) {
        vb[0] = vld1q_u16(&b[j * PARAMS_NBAR]);
        vb[1] = vld1q_u16(&b[(j + 1) * PARAMS_NBAR]);
        vb[2] = vld1q_u16(&b[(j + 2) * PARAMS_NBAR]);
        vb[3] = vld1q_u16(&b[(j + 3) * PARAMS_NBAR]);
        vb[4] = vld1q_u16(&b[(j + 4) * PARAMS_NBAR]);
        vb[5] = vld1q_u16(&b[(j + 5) * PARAMS_NBAR]);
        vb[6] = vld1q_u16(&b[(j + 6) * PARAMS_NBAR]);
        vb[7] = vld1q_u16(&b[(j + 7) * PARAMS_NBAR]);

        for (i = 0; i < PARAMS_NBAR; i++) {
            vs = vld1q_u16(&s[i * PARAMS_N + j]);

            vout[i] = vmlaq_laneq_u16(vout[i], vb[0], vs, 0);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[1], vs, 1);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[2], vs, 2);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[3], vs, 3);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[4], vs, 4);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[5], vs, 5);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[6], vs, 6);
            vout[i] = vmlaq_laneq_u16(vout[i], vb[7], vs, 7);
        }
    }

    for (i = 0; i < 8; i++) {
        vout[i] = vandq_u16(vout[i], vmask);
        vst1q_u16(&out[i * PARAMS_NBAR], vout[i]);
    }
#else
    int k;

    for (k = 0; k < PARAMS_NBAR; k++) {
        for (i = 0; i < PARAMS_NBAR; i++) {
            out[k*PARAMS_NBAR + i] = e[k*PARAMS_NBAR + i];
            for (j = 0; j < PARAMS_N; j++) {
                out[k*PARAMS_NBAR + i] += (int16_t)s[k*PARAMS_N + j] * b[j*PARAMS_NBAR + i];
            }
            out[k*PARAMS_NBAR + i] = (uint32_t)(out[k*PARAMS_NBAR + i]) & ((1<<PARAMS_LOGQ)-1);
        }
    }
#endif
}


void frodo_add(uint16_t *out, const uint16_t *a, const uint16_t *b) 
{ // Add a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a + b

    for (int i = 0; i < (PARAMS_NBAR*PARAMS_NBAR); i++) {
        out[i] = (a[i] + b[i]) & ((1<<PARAMS_LOGQ)-1);
    }
}


void frodo_sub(uint16_t *out, const uint16_t *a, const uint16_t *b) 
{ // Subtract a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a - b

    for (int i = 0; i < (PARAMS_NBAR*PARAMS_NBAR); i++) {
        out[i] = (a[i] - b[i]) & ((1<<PARAMS_LOGQ)-1);
    }
}


void frodo_key_encode(uint16_t *out, const uint16_t *in) 
{ // Encoding
    unsigned int i, j, npieces_word = 8;
    unsigned int nwords = (PARAMS_NBAR*PARAMS_NBAR)/8;
    uint64_t temp, mask = ((uint64_t)1 << PARAMS_EXTRACTED_BITS) - 1;
    uint16_t* pos = out;

    for (i = 0; i < nwords; i++) {
        temp = 0;
        for(j = 0; j < PARAMS_EXTRACTED_BITS; j++) 
            temp |= ((uint64_t)((uint8_t*)in)[i*PARAMS_EXTRACTED_BITS + j]) << (8*j);
        for (j = 0; j < npieces_word; j++) { 
            *pos = (uint16_t)((temp & mask) << (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS));  
            temp >>= PARAMS_EXTRACTED_BITS;
            pos++;
        }
    }
}


void frodo_key_decode(uint16_t *out, const uint16_t *in)
{ // Decoding
    unsigned int i, j, index = 0, npieces_word = 8;
    unsigned int nwords = (PARAMS_NBAR * PARAMS_NBAR) / 8;
    uint16_t temp, maskex=((uint16_t)1 << PARAMS_EXTRACTED_BITS) -1, maskq =((uint16_t)1 << PARAMS_LOGQ) -1;
    uint8_t  *pos = (uint8_t*)out;
    uint64_t templong;

    for (i = 0; i < nwords; i++) {
        templong = 0;
        for (j = 0; j < npieces_word; j++) {  // temp = floor(in*2^{-11}+0.5)
            temp = ((in[index] & maskq) + (1 << (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS - 1))) >> (PARAMS_LOGQ - PARAMS_EXTRACTED_BITS);
            templong |= ((uint64_t)(temp & maskex)) << (PARAMS_EXTRACTED_BITS * j);
            index++;
        }
	for(j = 0; j < PARAMS_EXTRACTED_BITS; j++) 
	    pos[i*PARAMS_EXTRACTED_BITS + j] = (templong >> (8*j)) & 0xFF;
    }
}
