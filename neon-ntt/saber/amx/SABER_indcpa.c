
/*
 * CC0 1.0 Universal or the following MIT License
 *
 * MIT License
 *
 * Copyright (c) 2023: Hanno Becker, Vincent Hwang, Matthias J. Kannwischer, Bo-Yin Yang, and Shang-Yi Yang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "SABER_indcpa.h"

#include <arm_neon.h>
#include <string.h>

#include "SABER_params.h"
#include "cbd.h"
#include "fips202.h"
#include "fips202x2.h"
#include "pack_unpack.h"
#include "poly.h"
#include "polymodmul.h"
#include "randombytes.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))

#ifdef ALLOC_MMAP
static uint16_t *A_kp, *s_kp, *b_kp;

__attribute__((constructor)) static void alloc_kp(void) {
    A_kp = MEMORY_ALLOC(SABER_L * SABER_L * SABER_N * sizeof(uint16_t));

    s_kp = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    b_kp = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
}
#endif

void indcpa_kem_keypair(uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES]) {
#ifdef ALLOC_STACK
    uint16_t A[SABER_L][SABER_L][SABER_N];

    uint16_t s[SABER_L][SABER_N];
    uint16_t b[SABER_L][SABER_N] = {0};
#else
    uint16_t(*A)[SABER_L][SABER_N] = (uint16_t(*)[SABER_L][SABER_N])A_kp;

    uint16_t(*s)[SABER_N] = (uint16_t(*)[SABER_N])s_kp;
    uint16_t(*b)[SABER_N] = (uint16_t(*)[SABER_N])b_kp;

    memset(b, 0, SABER_L * SABER_N * sizeof(uint16_t));
#endif

    uint8_t seed_A[SABER_SEEDBYTES];
    uint8_t seed_s[SABER_NOISE_SEEDBYTES];

    uint8_t shake_A_buf[SABER_L * SABER_L * SABER_POLYBYTES];
    uint8_t shake_s_buf[SABER_L * SABER_POLYCOINBYTES];

    randombytes(seed_A, SABER_SEEDBYTES);
    shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);  // for not revealing system RNG state
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);

#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || __ARM_FEATURE_SHA3
    uint8_t shake_A_buf_dummy[SABER_L * SABER_POLYVECBYTES];
    uint8_t shake_s_buf_dummy[SABER_L * SABER_POLYCOINBYTES];
    shake128x2(shake_A_buf, shake_A_buf_dummy, sizeof(shake_A_buf), seed_A, seed_A, SABER_SEEDBYTES);
    shake128x2(shake_s_buf, shake_s_buf_dummy, sizeof(shake_s_buf), seed_s, seed_s, SABER_NOISE_SEEDBYTES);
#else
    shake128(shake_A_buf, sizeof(shake_A_buf), seed_A, SABER_SEEDBYTES);
    shake128(shake_s_buf, sizeof(shake_s_buf), seed_s, SABER_NOISE_SEEDBYTES);
#endif

    for (int i = 0; i < SABER_L; i++) {
        BS2POLVECq(shake_A_buf + i * SABER_POLYVECBYTES, A[i]);
    }

    for (int i = 0; i < SABER_L; i++) {
        cbd(s[i], shake_s_buf + i * SABER_POLYCOINBYTES);
    }

    MatrixVectorMulRound((const uint16_t(*)[SABER_L][SABER_N])A, (const uint16_t(*)[SABER_N])s, (uint16_t(*)[SABER_N])b,
                         1);

    POLVECq2BS(sk, (const uint16_t(*)[SABER_N])s);
    POLVECp2BS(pk, (const uint16_t(*)[SABER_N])b);
    memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed_A, sizeof(seed_A));
}

#ifdef ALLOC_MMAP
static uint16_t *A_enc, *sp_enc, *bp_enc, *vp_enc, *b_enc;

__attribute__((constructor)) static void alloc_enc(void) {
    A_enc = MEMORY_ALLOC(SABER_L * SABER_L * SABER_N * sizeof(uint16_t));

    sp_enc = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    bp_enc = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    vp_enc = MEMORY_ALLOC(SABER_N * sizeof(uint16_t));
    b_enc = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
}
#endif

void indcpa_kem_enc(const uint8_t m[SABER_KEYBYTES], const uint8_t seed_sp[SABER_NOISE_SEEDBYTES],
                    const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t ciphertext[SABER_BYTES_CCA_DEC]) {
#ifdef ALLOC_STACK
    uint16_t A[SABER_L][SABER_L][SABER_N];

    uint16_t sp[SABER_L][SABER_N];
    uint16_t bp[SABER_L][SABER_N] = {0};
    uint16_t vp[SABER_N] = {0};
    uint16_t b[SABER_L][SABER_N];
#else
    uint16_t(*A)[SABER_L][SABER_N] = (uint16_t(*)[SABER_L][SABER_N])A_enc;

    uint16_t(*sp)[SABER_N] = (uint16_t(*)[SABER_N])sp_enc;
    uint16_t(*bp)[SABER_N] = (uint16_t(*)[SABER_N])bp_enc;
    uint16_t(*vp) = (uint16_t(*))vp_enc;
    uint16_t(*b)[SABER_N] = (uint16_t(*)[SABER_N])b_enc;
#endif

    uint16_t mp[SABER_N];
    const uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;

    uint8_t shake_A_buf[SABER_L * SABER_L * SABER_POLYBYTES];
    uint8_t shake_s_buf[SABER_L * SABER_POLYCOINBYTES];

#if (__APPLE__ && __ARM_FEATURE_CRYPTO) || __ARM_FEATURE_SHA3
    uint8_t shake_A_buf_dummy[SABER_L * SABER_POLYVECBYTES];
    uint8_t shake_s_buf_dummy[SABER_L * SABER_POLYCOINBYTES];
    shake128x2(shake_A_buf, shake_A_buf_dummy, sizeof(shake_A_buf), seed_A, seed_A, SABER_SEEDBYTES);
    shake128x2(shake_s_buf, shake_s_buf_dummy, sizeof(shake_s_buf), seed_sp, seed_sp, SABER_NOISE_SEEDBYTES);
#else
    shake128(shake_A_buf, sizeof(shake_A_buf), seed_A, SABER_SEEDBYTES);
    shake128(shake_s_buf, sizeof(shake_s_buf), seed_sp, SABER_NOISE_SEEDBYTES);
#endif

    for (int i = 0; i < SABER_L; i++) {
        BS2POLVECq(shake_A_buf + i * SABER_POLYVECBYTES, A[i]);
    }

    for (int i = 0; i < SABER_L; i++) {
        cbd(sp[i], shake_s_buf + i * SABER_POLYCOINBYTES);
    }

    MatrixVectorMulRound((const uint16_t(*)[SABER_L][SABER_N])A, (const uint16_t(*)[SABER_N])sp,
                         (uint16_t(*)[SABER_N])bp, 0);

    BS2POLVECp(pk, b);
    BS2POLmsg(m, mp);

    InnerProd((const uint16_t(*)[SABER_N])b, (const uint16_t(*)[SABER_N])sp, vp);

    for (int j = 0; j < SABER_N; j++) {
        vp[j] = (vp[j] - (mp[j] << (SABER_EP - 1)) + h1) >> (SABER_EP - SABER_ET);
    }

    POLVECp2BS(ciphertext, (const uint16_t(*)[SABER_N])bp);
    POLT2BS(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, vp);
}

#ifdef ALLOC_MMAP
static uint16_t *b_dec, *s_dec, *v_dec;

__attribute__((constructor)) static void alloc_dec(void) {
    b_dec = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));

    s_dec = MEMORY_ALLOC(SABER_L * SABER_N * sizeof(uint16_t));
    v_dec = MEMORY_ALLOC(SABER_N * sizeof(uint16_t));
}
#endif

void indcpa_kem_dec(const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC],
                    uint8_t m[SABER_KEYBYTES]) {
#ifdef ALLOC_STACK
    uint16_t b[SABER_L][SABER_N];
    uint16_t s[SABER_L][SABER_N];

    uint16_t v[SABER_N] = {0};
#else
    uint16_t(*b)[SABER_N] = (uint16_t(*)[SABER_N])b_dec;
    uint16_t(*s)[SABER_N] = (uint16_t(*)[SABER_N])s_dec;

    uint16_t *v = v_dec;
#endif
    uint16_t cm[SABER_N];

    BS2POLT(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, cm);

    BS2POLVECq(sk, s);
    BS2POLVECp(ciphertext, b);

    InnerProd((const uint16_t(*)[SABER_N])b, (const uint16_t(*)[SABER_N])s, v);

    for (int i = 0; i < SABER_N; i++) {
        v[i] = (v[i] + h2 - (cm[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
    }

    POLmsg2BS(m, v);
}
