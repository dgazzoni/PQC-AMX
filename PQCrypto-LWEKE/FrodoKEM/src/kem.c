/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: Key Encapsulation Mechanism (KEM) based on Frodo
*********************************************************************************************/

#include <string.h>
#include "fips202.h"
#include "../../common/random/random.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#ifdef DO_VALGRIND_CHECK
#include <valgrind/memcheck.h>
#endif

#ifdef ALLOC_MMAP
uint16_t *B_kp, *S_kp;

__attribute__((constructor)) static void alloc_kp(void) {
    B_kp = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    S_kp = MEMORY_ALLOC(2*PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
}
#endif

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk)
{ // FrodoKEM's key generation
  // Outputs: public key pk = pk_seedA||pk_b                      (               BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 bytes)
  //          secret key sk = sk_s||pk_seedA||pk_b||sk_S||sk_pkh  (CRYPTO_BYTES + BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + 2*PARAMS_N*PARAMS_NBAR + BYTES_PKHASH bytes)
    uint8_t *pk_seedA = &pk[0];
    uint8_t *pk_b = &pk[BYTES_SEED_A];
    uint8_t *sk_s = &sk[0];
    uint8_t *sk_pk = &sk[CRYPTO_BYTES];
    uint8_t *sk_S = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES];
    uint8_t *sk_pkh = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES + 2*PARAMS_N*PARAMS_NBAR];
#ifdef ALLOC_MMAP
    uint16_t *B = B_kp;
    uint16_t *S = S_kp;
#else
    ALIGN_HEADER(128) uint16_t B[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};
    ALIGN_HEADER(128) uint16_t S[2*PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(128) = {0}; // contains secret data
#endif
    uint16_t *E = (uint16_t *)&S[PARAMS_N*PARAMS_NBAR];                // contains secret data
    uint8_t randomness[CRYPTO_BYTES + BYTES_SEED_SE + BYTES_SEED_A];   // contains secret data via randomness_s and randomness_seedSE
    uint8_t *randomness_s = &randomness[0];                            // contains secret data
    uint8_t *randomness_seedSE = &randomness[CRYPTO_BYTES];            // contains secret data
    uint8_t *randomness_z = &randomness[CRYPTO_BYTES + BYTES_SEED_SE];
    uint8_t shake_input_seedSE[1 + BYTES_SEED_SE];                     // contains secret data

#ifdef ALLOC_MMAP
    memset(B, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    memset(S, 0, 2*PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
#endif

    // Generate the secret value s, the seed for S and E, and the seed for the seed for A. Add seed_A to the public key
    if (randombytes(randomness, CRYPTO_BYTES + BYTES_SEED_SE + BYTES_SEED_A) != 0)
        return 1;
#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_UNDEFINED(randomness, CRYPTO_BYTES + BYTES_SEED_SE + BYTES_SEED_A);
#endif
    shake(pk_seedA, BYTES_SEED_A, randomness_z, BYTES_SEED_A);

    // Generate S and E, and compute B = A*S + E. Generate A on-the-fly
    shake_input_seedSE[0] = 0x5F;
    memcpy(&shake_input_seedSE[1], randomness_seedSE, BYTES_SEED_SE);
    shake((uint8_t*)S, 2*PARAMS_N*PARAMS_NBAR*sizeof(uint16_t), shake_input_seedSE, 1 + BYTES_SEED_SE);
    for (size_t i = 0; i < 2 * PARAMS_N * PARAMS_NBAR; i++) {
        S[i] = LE_TO_UINT16(S[i]);
    }
    frodo_sample_n(S, PARAMS_N*PARAMS_NBAR);
    frodo_sample_n(E, PARAMS_N*PARAMS_NBAR);
    frodo_mul_add_as_plus_e(B, S, E, pk);

    // Encode the second part of the public key
    frodo_pack(pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, B, PARAMS_N*PARAMS_NBAR, PARAMS_LOGQ);

    // Add s, pk and S to the secret key
    memcpy(sk_s, randomness_s, CRYPTO_BYTES);
    memcpy(sk_pk, pk, CRYPTO_PUBLICKEYBYTES);
    for (size_t i = 0; i < PARAMS_N * PARAMS_NBAR; i++) {
        S[i] = UINT16_TO_LE(S[i]);
    }
    memcpy(sk_S, S, 2*PARAMS_N*PARAMS_NBAR);

    // Add H(pk) to the secret key
    shake(sk_pkh, BYTES_PKHASH, pk, CRYPTO_PUBLICKEYBYTES);

    // Cleanup:
    clear_bytes((uint8_t *)S, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)E, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes(randomness, CRYPTO_BYTES + BYTES_SEED_SE);
    clear_bytes(shake_input_seedSE, 1 + BYTES_SEED_SE);
#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_DEFINED(randomness, CRYPTO_BYTES + BYTES_SEED_SE + BYTES_SEED_A);
#endif
    return 0;
}

#ifdef ALLOC_MMAP
uint16_t *Bp_enc, *Sp_enc, *Ep_enc, *Epp_enc;

__attribute__((constructor)) static void alloc_enc(void) {
    Bp_enc = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    Sp_enc = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    Ep_enc = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    Epp_enc = MEMORY_ALLOC(PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
}

#define CONCAT2(x, y) x##y
#define CONCAT(x, y) CONCAT2(x, y)
#define shakeincctx CONCAT(shake, incctx)
#define shake_inc_init CONCAT(shake, _inc_init)
#define shake_inc_absorb CONCAT(shake, _inc_absorb)
#define shake_inc_finalize CONCAT(shake, _inc_finalize)
#define shake_inc_squeeze CONCAT(shake, _inc_squeeze)
#define shake_inc_ctx_release CONCAT(shake, _inc_ctx_release)
#endif

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{ // FrodoKEM's key encapsulation
  // Input:   public key pk = pk_seedA||pk_b      (BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 bytes)
  // Outputs: ciphertext ct = ct_c1||ct_c2||salt  (               (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8 + BYTES_SALT bytes)
  //          shared key ss                       (CRYPTO_BYTES bytes)
    const uint8_t *pk_seedA = &pk[0];
    const uint8_t *pk_b = &pk[BYTES_SEED_A];
    uint8_t *ct_c1 = &ct[0];
    uint8_t *ct_c2 = &ct[(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8];
    uint16_t B[PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t V[PARAMS_NBAR*PARAMS_NBAR]= {0};                          // contains secret data
    uint16_t C[PARAMS_NBAR*PARAMS_NBAR] = {0};
#ifdef ALLOC_MMAP
    uint16_t *Bp = Bp_enc;
    // The combined allocation of Sp, Ep and Epp, as done in the original code, leads to significant performance
    // degradation for mmap() allocation. This was fixed through separate allocations for each.
    uint16_t *Sp = Sp_enc;                                             // contains secret data
    uint16_t *Ep = Ep_enc;                                             // contains secret data
    uint16_t *Epp = Epp_enc;                                           // contains secret data
#else
    ALIGN_HEADER(32) uint16_t Bp[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(128) uint16_t Sp[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data
    uint16_t *Ep = (uint16_t *)&Sp[PARAMS_N*PARAMS_NBAR];              // contains secret data
    uint16_t *Epp = (uint16_t *)&Sp[2*PARAMS_N*PARAMS_NBAR];           // contains secret data
#endif
    uint8_t G2in[BYTES_PKHASH + BYTES_MU + BYTES_SALT];                // contains secret data via mu
    uint8_t *pkh = &G2in[0];
    uint8_t *mu = &G2in[BYTES_PKHASH];                                 // contains secret data
    uint8_t *salt = &G2in[BYTES_PKHASH + BYTES_MU];
    uint8_t G2out[BYTES_SEED_SE + CRYPTO_BYTES];                       // contains secret data
    uint8_t *seedSE = &G2out[0];                                       // contains secret data
    uint8_t *k = &G2out[BYTES_SEED_SE];                                // contains secret data
    uint8_t Fin[CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];                // contains secret data via Fin_k
    uint8_t *Fin_ct = &Fin[0];
    uint8_t *Fin_k = &Fin[CRYPTO_CIPHERTEXTBYTES];                     // contains secret data
    uint8_t shake_input_seedSE[1 + BYTES_SEED_SE];                     // contains secret data

#ifdef ALLOC_MMAP
    memset(Bp, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    memset(Sp, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    memset(Ep, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    memset(Epp, 0, PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
#endif

    // pkh <- G_1(pk), generate random mu and salt, compute (seedSE || k) = G_2(pkh || mu || salt)
    shake(pkh, BYTES_PKHASH, pk, CRYPTO_PUBLICKEYBYTES);
    if (randombytes(mu, BYTES_MU + BYTES_SALT) != 0)
        return 1;
#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_UNDEFINED(mu, BYTES_MU + BYTES_SALT);
    VALGRIND_MAKE_MEM_UNDEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif
    shake(G2out, BYTES_SEED_SE + CRYPTO_BYTES, G2in, BYTES_PKHASH + BYTES_MU + BYTES_SALT);

    // Generate Sp and Ep, and compute Bp = Sp*A + Ep. Generate A on-the-fly
    shake_input_seedSE[0] = 0x96;
    memcpy(&shake_input_seedSE[1], seedSE, BYTES_SEED_SE);

#ifdef ALLOC_MMAP
    shakeincctx ctx;

    shake_inc_init(&ctx);
    shake_inc_absorb(&ctx, shake_input_seedSE, 1 + BYTES_SEED_SE);
    shake_inc_finalize(&ctx);
    shake_inc_squeeze((uint8_t*)Sp, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t), &ctx);
    shake_inc_squeeze((uint8_t*)Ep, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t), &ctx);
    shake_inc_squeeze((uint8_t*)Epp, PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t), &ctx);
    shake_inc_ctx_release(&ctx);
#else
    shake((uint8_t*)Sp, (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t), shake_input_seedSE, 1 + BYTES_SEED_SE);
#endif
    frodo_sample_n(Sp, PARAMS_N*PARAMS_NBAR);
    frodo_sample_n(Ep, PARAMS_N*PARAMS_NBAR);
    frodo_mul_add_sa_plus_e(Bp, Sp, Ep, pk_seedA);
    frodo_pack(ct_c1, (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8, Bp, PARAMS_N*PARAMS_NBAR, PARAMS_LOGQ);

    // Generate Epp, and compute V = Sp*B + Epp
    frodo_sample_n(Epp, PARAMS_NBAR*PARAMS_NBAR);
    frodo_unpack(B, PARAMS_N*PARAMS_NBAR, pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, PARAMS_LOGQ);
    frodo_mul_add_sb_plus_e(V, B, Sp, Epp);

    // Encode mu, and compute C = V + enc(mu) (mod q)
    frodo_key_encode(C, (uint16_t*)mu);
    frodo_add(C, V, C);
    frodo_pack(ct_c2, (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8, C, PARAMS_NBAR*PARAMS_NBAR, PARAMS_LOGQ);

    // Append salt to ct and compute ss = F(ct_c1||ct_c2||salt||k)
    memcpy(&ct[CRYPTO_CIPHERTEXTBYTES - BYTES_SALT], salt, BYTES_SALT);
    memcpy(Fin_ct, ct, CRYPTO_CIPHERTEXTBYTES);
    memcpy(Fin_k, k, CRYPTO_BYTES);
    shake(ss, CRYPTO_BYTES, Fin, CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);

    // Cleanup:
    clear_bytes((uint8_t *)V, PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Sp, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Ep, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Epp, PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes(mu, BYTES_MU);
    clear_bytes(G2out, BYTES_SEED_SE + CRYPTO_BYTES);
    clear_bytes(Fin_k, CRYPTO_BYTES);
    clear_bytes(shake_input_seedSE, 1 + BYTES_SEED_SE);
#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_DEFINED(mu, BYTES_MU);
    VALGRIND_MAKE_MEM_DEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif
    return 0;
}

#ifdef ALLOC_MMAP
uint16_t *BBp_dec, *Sp_dec;

__attribute__((constructor)) static void alloc_dec(void) {
    BBp_dec = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    Sp_dec = MEMORY_ALLOC((2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR * sizeof(uint16_t));
}
#endif

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{ // FrodoKEM's key decapsulation
  // Inputs: ciphertext ct = ct_c1||ct_c2||salt                  (                              (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8 + BYTES_SALT bytes)
  //         secret key sk = sk_s||pk_seedA||pk_b||sk_S||sk_pkh  (CRYPTO_BYTES + BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + 2*PARAMS_N*PARAMS_NBAR + BYTES_PKHASH bytes)
  // Output: shared key ss                                       (CRYPTO_BYTES bytes)
    uint16_t B[PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t Bp[PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t W[PARAMS_NBAR*PARAMS_NBAR] = {0};                         // contains secret data
    uint16_t C[PARAMS_NBAR*PARAMS_NBAR] = {0};
    uint16_t CC[PARAMS_NBAR*PARAMS_NBAR] = {0};
#ifdef ALLOC_MMAP
    uint16_t *BBp = BBp_dec;
    uint16_t *Sp = Sp_dec;
#else
    ALIGN_HEADER(32) uint16_t BBp[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(128) uint16_t Sp[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data
#endif
    uint16_t *Ep = (uint16_t *)&Sp[PARAMS_N*PARAMS_NBAR];              // contains secret data
    uint16_t *Epp = (uint16_t *)&Sp[2*PARAMS_N*PARAMS_NBAR];           // contains secret data
    const uint8_t *ct_c1 = &ct[0];
    const uint8_t *ct_c2 = &ct[(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8];
    const uint8_t *salt = &ct[CRYPTO_CIPHERTEXTBYTES - BYTES_SALT];
    const uint8_t *sk_s = &sk[0];
    const uint8_t *sk_pk = &sk[CRYPTO_BYTES];
    const uint16_t *sk_S = (uint16_t *) &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES];
    uint16_t S[PARAMS_N * PARAMS_NBAR];                                // contains secret data
    const uint8_t *sk_pkh = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES + 2*PARAMS_N*PARAMS_NBAR];
    const uint8_t *pk_seedA = &sk_pk[0];
    const uint8_t *pk_b = &sk_pk[BYTES_SEED_A];
    uint8_t G2in[BYTES_PKHASH + BYTES_MU + BYTES_SALT];                // contains secret data via muprime
    uint8_t *pkh = &G2in[0];
    uint8_t *muprime = &G2in[BYTES_PKHASH];                            // contains secret data
    uint8_t *G2in_salt = &G2in[BYTES_PKHASH + BYTES_MU];
    uint8_t G2out[BYTES_SEED_SE + CRYPTO_BYTES];                       // contains secret data
    uint8_t *seedSEprime = &G2out[0];                                  // contains secret data
    uint8_t *kprime = &G2out[BYTES_SEED_SE];                           // contains secret data
    uint8_t Fin[CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];                // contains secret data via Fin_k
    uint8_t *Fin_ct = &Fin[0];
    uint8_t *Fin_k = &Fin[CRYPTO_CIPHERTEXTBYTES];                     // contains secret data
    uint8_t shake_input_seedSEprime[1 + BYTES_SEED_SE];                // contains secret data

#ifdef ALLOC_MMAP
    memset(BBp, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    memset(Sp, 0, (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR * sizeof(uint16_t));
#endif

#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_UNDEFINED(sk, CRYPTO_SECRETKEYBYTES);
    VALGRIND_MAKE_MEM_UNDEFINED(ct, CRYPTO_CIPHERTEXTBYTES);
#endif

    for (size_t i = 0; i < PARAMS_N * PARAMS_NBAR; i++) {
        S[i] = LE_TO_UINT16(sk_S[i]);
    }

    // Compute W = C - Bp*S (mod q), and decode the randomness mu
    frodo_unpack(Bp, PARAMS_N*PARAMS_NBAR, ct_c1, (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8, PARAMS_LOGQ);
    frodo_unpack(C, PARAMS_NBAR*PARAMS_NBAR, ct_c2, (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8, PARAMS_LOGQ);
    frodo_mul_bs(W, Bp, S);
    frodo_sub(W, C, W);
    frodo_key_decode((uint16_t*)muprime, W);

    // Generate (seedSE' || k') = G_2(pkh || mu' || salt)
    memcpy(pkh, sk_pkh, BYTES_PKHASH);
    memcpy(G2in_salt, salt, BYTES_SALT);
    shake(G2out, BYTES_SEED_SE + CRYPTO_BYTES, G2in, BYTES_PKHASH + BYTES_MU + BYTES_SALT);

    // Generate Sp and Ep, and compute BBp = Sp*A + Ep. Generate A on-the-fly
    shake_input_seedSEprime[0] = 0x96;
    memcpy(&shake_input_seedSEprime[1], seedSEprime, BYTES_SEED_SE);
    shake((uint8_t*)Sp, (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t), shake_input_seedSEprime, 1 + BYTES_SEED_SE);
    for (size_t i = 0; i < (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR; i++) {
        Sp[i] = LE_TO_UINT16(Sp[i]);
    }
    frodo_sample_n(Sp, PARAMS_N*PARAMS_NBAR);
    frodo_sample_n(Ep, PARAMS_N*PARAMS_NBAR);
    frodo_mul_add_sa_plus_e(BBp, Sp, Ep, pk_seedA);

    // Generate Epp, and compute W = Sp*B + Epp
    frodo_sample_n(Epp, PARAMS_NBAR*PARAMS_NBAR);
    frodo_unpack(B, PARAMS_N*PARAMS_NBAR, pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, PARAMS_LOGQ);
    frodo_mul_add_sb_plus_e(W, B, Sp, Epp);

    // Encode mu, and compute CC = W + enc(mu') (mod q)
    frodo_key_encode(CC, (uint16_t*)muprime);
    frodo_add(CC, W, CC);

    // Prepare input to F
    memcpy(Fin_ct, ct, CRYPTO_CIPHERTEXTBYTES);

    // Reducing BBp modulo q
    for (int i = 0; i < PARAMS_N*PARAMS_NBAR; i++) BBp[i] = BBp[i] & ((1 << PARAMS_LOGQ)-1);

    // If (Bp == BBp & C == CC) then ss = F(ct || k'), else ss = F(ct || s)
    // Needs to avoid branching on secret data using constant-time implementation.
    int8_t selector = ct_verify(Bp, BBp, PARAMS_N*PARAMS_NBAR) | ct_verify(C, CC, PARAMS_NBAR*PARAMS_NBAR);
    // If (selector == 0) then load k' to do ss = F(ct || k'), else if (selector == -1) load s to do ss = F(ct || s)
    ct_select((uint8_t*)Fin_k, (uint8_t*)kprime, (uint8_t*)sk_s, CRYPTO_BYTES, selector);
    shake(ss, CRYPTO_BYTES, Fin, CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);

    // Cleanup:
    clear_bytes((uint8_t *)W, PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Sp, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)S, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Ep, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes((uint8_t *)Epp, PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
    clear_bytes(muprime, BYTES_MU);
    clear_bytes(G2out, BYTES_SEED_SE + CRYPTO_BYTES);
    clear_bytes(Fin_k, CRYPTO_BYTES);
    clear_bytes(shake_input_seedSEprime, 1 + BYTES_SEED_SE);
#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_DEFINED(sk, CRYPTO_SECRETKEYBYTES);
    VALGRIND_MAKE_MEM_DEFINED(ct, CRYPTO_CIPHERTEXTBYTES);
#endif
    return 0;
}
