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
#include "frodokem_opt_matmul.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

void AES128_ECB_genA_rows(const uint8_t *schedule, uint16_t *A, int n_rows, int first_row, int n, int n_cols);
void AES128_ECB_genA_32cols(const uint8_t *schedule, uint16_t *A, int n, int first_col);

#ifdef ALLOC_MMAP
int16_t *a_row_as;
uint16_t *s_as;

__attribute__((constructor)) static void alloc_as(void) {
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    a_row_as = MEMORY_ALLOC(32 * n32 * sizeof(int16_t));
    s_as = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_as_plus_e(uint16_t *out, const uint16_t *sT, const uint16_t *e, const uint8_t *seed_A) 
{ // Generate-and-multiply: generate matrix A (N x N) row-wise, multiply by s on the right.
  // Inputs: s, e (N x N_BAR)
  // Output: out = A*s + e (N x N_BAR)
    int i, j;
    const int n32 = 32 * ((PARAMS_N + 31) / 32);

#ifdef ALLOC_MMAP
    int16_t *a_row = a_row_as;
    uint16_t *s = s_as;
#else
    ALIGN_HEADER(32) int16_t a_row[32 * n32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif
    memset(a_row, 0, 32 * n32 * sizeof(int16_t));
    
    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
        }
    }

#if defined(USE_AES128_FOR_A)
    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);   
                                    
    for (i = 0; i < PARAMS_N; i += 32) {
        AES128_ECB_genA_rows(aes_key_schedule, (uint16_t*)a_row, MIN(PARAMS_N - i, 32), i, PARAMS_N, n32);
#elif defined (USE_SHAKE128_FOR_A)       
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
            shake128x2((uint8_t*)(a_row + j*n32), (uint8_t*)(a_row + (j + 1)*n32), (unsigned long long)(2*PARAMS_N), seed_A_separated_0, seed_A_separated_1, 2 + BYTES_SEED_A);
        }
#endif

        amx_mul_add_as_plus_e_up_to_32_rows(out, a_row, s, e, PARAMS_N, i);
    }
    
#if defined(USE_AES128_FOR_A)
    AES128_free_schedule(aes_key_schedule);
#endif
    return 1;
}

#if defined(USE_AES128_FOR_A)
#ifdef ALLOC_MMAP
uint16_t *A_sa, *s_sa;

__attribute__((constructor)) static void alloc_sa(void) {
    A_sa = MEMORY_ALLOC(PARAMS_N * 32 * sizeof(int16_t));
    s_sa = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e(uint16_t *out, const uint16_t *sT, uint16_t *e, const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j;
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa;
    uint16_t *s = s_sa;
#else
    ALIGN_HEADER(32) uint16_t A[PARAMS_N * 32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif
    memset(A, 0, PARAMS_N * 32 * sizeof(uint16_t));

    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
        }
    }

    uint8_t aes_key_schedule[16*11];
    AES128_load_schedule(seed_A, aes_key_schedule);

    // Start matrix multiplication
    for (i = 0; i < PARAMS_N; i += 32) {
        AES128_ECB_genA_32cols(aes_key_schedule, A, i, PARAMS_N);
        amx_mul_add_sa_plus_e_up_to_32_cols(out, (int16_t *)A, s, e, PARAMS_N, i);
    }

    AES128_free_schedule(aes_key_schedule);

    return 1;
}
#elif defined (USE_SHAKE128_FOR_A)  // SHAKE128
#ifdef ALLOC_MMAP
uint16_t *A_sa, *s_sa;

__attribute__((constructor)) static void alloc_sa(void) {
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
    A_sa = MEMORY_ALLOC(32 * n32 * sizeof(int16_t));
    s_sa = MEMORY_ALLOC((PARAMS_N / 4) * 32 * sizeof(uint16_t));
}
#endif

int frodo_mul_add_sa_plus_e(uint16_t *out, const uint16_t *sT, uint16_t *e, const uint8_t *seed_A)
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
  // The matrix multiplication uses the row-wise blocking and packing (RWCF) approach described in: J.W. Bos, M. Ofner, J. Renes, 
  // T. Schneider, C. van Vredendaal, "The Matrix Reloaded: Multiplication Strategies in FrodoKEM". https://eprint.iacr.org/2021/711
    int i, j;
    const int n32 = 32 * ((PARAMS_N + 31) / 32);
#ifdef ALLOC_MMAP
    uint16_t *A = A_sa;
    uint16_t *s = s_sa;
#else
    ALIGN_HEADER(32) uint16_t A[32 * n32] ALIGN_FOOTER(32);
    ALIGN_HEADER(32) uint16_t s[(PARAMS_N / 4) * 32] ALIGN_FOOTER(32);
#endif
    memset(A, 0, 32 * n32 * sizeof(uint16_t));

    for (i = 0; i < PARAMS_N; i++) {
        for (j = 0; j < 8; j++) {
            s[i * 8 + j] = sT[j * PARAMS_N + i];
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

        amx_mul_add_sa_plus_e_up_to_32_rows(out, (int16_t *)A, s, e, PARAMS_N, i);
    }

    return 1;
}
#endif

void frodo_mul_bs(uint16_t *out, const uint16_t *b, const uint16_t *s) 
{ // Multiply by s on the right
  // Inputs: b (N_BAR x N), s (N x N_BAR)
  // Output: out = b*s (N_BAR x N_BAR)
    int i, j, k;

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
}


void frodo_mul_add_sb_plus_e(uint16_t *out, const uint16_t *b, const uint16_t *s, const uint16_t *e) 
{ // Multiply by s on the left
  // Inputs: b (N x N_BAR), s (N_BAR x N), e (N_BAR x N_BAR)
  // Output: out = s*b + e (N_BAR x N_BAR)
    int i, j;

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
