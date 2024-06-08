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
#else
    #include "fips202.h"
#endif
#endif    

#ifdef ALLOC_MMAP
int16_t *A_sa_4x;

__attribute__((constructor)) static void alloc_sa_4x(void) {
    A_sa_4x = MEMORY_ALLOC(PARAMS_N * PARAMS_N * sizeof(int16_t));
}
#endif

int frodo_mul_add_sa_plus_e_4x(uint16_t *out[4], const uint16_t *s[4], uint16_t *e[4], const uint8_t *seed_A) 
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
    int i, j, k, n;
#ifdef ALLOC_MMAP
    int16_t *A = A_sa_4x;
#else
    int16_t A[PARAMS_N * PARAMS_N];
#endif
    memset(A, 0, PARAMS_N * PARAMS_N * sizeof(int16_t));
    
#if defined(USE_AES128_FOR_A)    // Matrix A generation using AES128, done per 128-bit block                                       
    size_t A_len = PARAMS_N * PARAMS_N * sizeof(int16_t);      
    for (i = 0; i < PARAMS_N; i++) {                        
        for (j = 0; j < PARAMS_N; j += PARAMS_STRIPE_STEP) {
            A[i*PARAMS_N + j] = UINT16_TO_LE(i);                // Loading values in the little-endian order
            A[i*PARAMS_N + j + 1] = UINT16_TO_LE(j);
        }
    }
    
#if !defined(USE_OPENSSL)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);  
    AES128_ECB_enc_sch((uint8_t*)A, A_len, aes_key_schedule, (uint8_t*)A);
#else
    EVP_CIPHER_CTX *aes_key_schedule;    
    int len;
    if (!(aes_key_schedule = EVP_CIPHER_CTX_new())) handleErrors();    
    if (1 != EVP_EncryptInit_ex(aes_key_schedule, EVP_aes_128_ecb(), NULL, seed_A, NULL)) handleErrors();    
    if (1 != EVP_EncryptUpdate(aes_key_schedule, (uint8_t*)A, &len, (uint8_t*)A, A_len)) handleErrors();
#endif
#elif defined (USE_SHAKE128_FOR_A)  // Matrix A generation using SHAKE128, done per 16*N-bit row
#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || (__ARM_FEATURE_SHA3)
    uint8_t seed_A_separated_0[2 + BYTES_SEED_A];
    uint8_t seed_A_separated_1[2 + BYTES_SEED_A];
    uint16_t *seed_A_origin_0 = (uint16_t*)&seed_A_separated_0;
    uint16_t *seed_A_origin_1 = (uint16_t*)&seed_A_separated_1;
    memcpy(&seed_A_separated_0[2], seed_A, BYTES_SEED_A);
    memcpy(&seed_A_separated_1[2], seed_A, BYTES_SEED_A);

    for (i = 0; i < PARAMS_N; i += 2) {
        seed_A_origin_0[0] = UINT16_TO_LE(i + 0);
        seed_A_origin_1[0] = UINT16_TO_LE(i + 1);
        shake128x2((uint8_t*)(A + i*PARAMS_N), (uint8_t*)(A + (i + 1)*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
    }
#else
    uint8_t seed_A_separated[2 + BYTES_SEED_A];
    uint16_t* seed_A_origin = (uint16_t*)&seed_A_separated;
    memcpy(&seed_A_separated[2], seed_A, BYTES_SEED_A);
    for (i = 0; i < PARAMS_N; i++) {
        seed_A_origin[0] = UINT16_TO_LE((uint16_t) i);
        shake128((unsigned char*)(A + i*PARAMS_N), (unsigned long long)(2*PARAMS_N), seed_A_separated, 2 + BYTES_SEED_A);
    }
#endif
#endif
    for (i = 0; i < PARAMS_N * PARAMS_N; i++) {
        A[i] = LE_TO_UINT16(A[i]);
    }

    for (n = 0; n < 4; n++) {
        memcpy(out[n], e[n], PARAMS_NBAR * PARAMS_N * sizeof(uint16_t));

        for (i = 0; i < PARAMS_N; i++) {                            // Matrix multiplication-addition A*s + e
            for (k = 0; k < PARAMS_NBAR; k++) {
                uint16_t sum = 0;
                for (j = 0; j < PARAMS_N; j++) {                                
                    sum += A[j*PARAMS_N + i] * s[n][k*PARAMS_N + j];  
                }
                out[n][k*PARAMS_N + i] += sum;                      // Adding e. No need to reduce modulo 2^15, extra bits are taken care of during packing later on.
            }
        }
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
