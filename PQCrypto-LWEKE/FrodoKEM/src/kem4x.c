/********************************************************************************************
 * FrodoKEM: Learning with Errors Key Encapsulation
 *
 * Abstract: Key Encapsulation Mechanism (KEM) based on Frodo
 *********************************************************************************************/

#include <string.h>

#include "../../common/random/random.h"
#include "fips202.h"
#include "fips202x2.h"
#ifdef ALLOC_MMAP
#include "memory_alloc.h"
#endif

#ifdef DO_VALGRIND_CHECK
#include <valgrind/memcheck.h>
#endif

#ifdef ALLOC_MMAP
uint16_t *B_enc4x, *V_enc4x[4], *Bp_enc4x[4], *Sp_enc4x[4], *Ep_enc4x[4], *Epp_enc4x[4];

__attribute__((constructor)) static void alloc_enc4x(void) {
    B_enc4x = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    for (int i = 0; i < 4; i++) {
        V_enc4x[i] = MEMORY_ALLOC(PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
        Bp_enc4x[i] = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        Sp_enc4x[i] = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t) + SHAKE_RATE);
        Ep_enc4x[i] = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t) + SHAKE_RATE);
        Epp_enc4x[i] = MEMORY_ALLOC(PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t) + SHAKE_RATE);
    }
}

#define CONCAT2(x, y) x##y
#define CONCAT(x, y) CONCAT2(x, y)
#define shakex2_absorb CONCAT(shakex2, _absorb)
#define shakex2_squeezeblocks CONCAT(shakex2, _squeezeblocks)
#endif

int crypto_kem_enc4x(unsigned char *ct0, unsigned char *ss0, unsigned char *ct1, unsigned char *ss1, unsigned char *ct2,
                     unsigned char *ss2, unsigned char *ct3, unsigned char *ss3, const unsigned char *pk)
{ // FrodoKEM's key encapsulation
  // Input:   public key pk = pk_seedA||pk_b      (BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 bytes)
  // Outputs: ciphertext ct = ct_c1||ct_c2||salt  (               (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8 + BYTES_SALT bytes)
  //          shared key ss                       (CRYPTO_BYTES bytes)
    unsigned char* ct[4] = {ct0, ct1, ct2, ct3};
    unsigned char* ss[4] = {ss0, ss1, ss2, ss3};
    const uint8_t *pk_seedA = &pk[0];
    const uint8_t *pk_b = &pk[BYTES_SEED_A];
    uint8_t *ct_c1[4] = { ct[0], ct[1], ct[2], ct[3] };
    uint8_t *ct_c2[4] = {
        &ct[0][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[1][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[2][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[3][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8]
    };
    uint16_t C[4][PARAMS_NBAR*PARAMS_NBAR] = {0};
#ifdef ALLOC_MMAP
    uint16_t *B = B_enc4x;

    uint16_t *V[4] = { V_enc4x[0], V_enc4x[1], V_enc4x[2], V_enc4x[3] }; // contains secret data

    uint16_t *Bp[4] = { Bp_enc4x[0], Bp_enc4x[1], Bp_enc4x[2], Bp_enc4x[3] };
    // The combined allocation of Sp, Ep and Epp, as done in the original code, leads to significant performance
    // degradation for mmap() allocation. This was fixed through separate allocations for each.
    uint16_t *Sp[4] = {                                                // contains secret data
        Sp_enc4x[0], Sp_enc4x[1], Sp_enc4x[2], Sp_enc4x[3]
    };
    uint16_t *Ep[4] = {                                                // contains secret data
        Ep_enc4x[0], Ep_enc4x[1], Ep_enc4x[2], Ep_enc4x[3]
    };
    uint16_t *Epp[4] = {                                               // contains secret data
        Epp_enc4x[0], Epp_enc4x[1], Epp_enc4x[2], Epp_enc4x[3]
    };
#else
    uint16_t B[PARAMS_N*PARAMS_NBAR] = {0};

    ALIGN_HEADER(32) uint16_t V0[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t V1[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t V2[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t V3[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};

    uint16_t *V[4] = { V0, V1, V2, V3 };

    ALIGN_HEADER(32) uint16_t Bp0[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t Bp1[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t Bp2[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t Bp3[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};

    uint16_t *Bp[4] = { Bp0, Bp1, Bp2, Bp3 };

    ALIGN_HEADER(128) uint16_t Sp0[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data
    ALIGN_HEADER(128) uint16_t Sp1[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data
    ALIGN_HEADER(128) uint16_t Sp2[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data
    ALIGN_HEADER(128) uint16_t Sp3[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};  // contains secret data

    uint16_t *Sp[4] = { Sp0, Sp1, Sp2, Sp3 };                          // contains secret data

    uint16_t *Ep[4] = {                                                // contains secret data
        (uint16_t *)&Sp[0][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[1][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[2][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[3][PARAMS_N*PARAMS_NBAR]
    };
    uint16_t *Epp[4] = {                                               // contains secret data
        (uint16_t *)&Sp[0][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[1][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[2][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[3][2*PARAMS_N*PARAMS_NBAR]
    };
#endif
    uint8_t G2in[4][BYTES_PKHASH + BYTES_MU + BYTES_SALT];             // contains secret data via mu
    uint8_t *pkh[4] = { G2in[0], G2in[1], G2in[2], G2in[3] };
    uint8_t *mu[4] = {                                                 // contains secret data
        &G2in[0][BYTES_PKHASH], &G2in[1][BYTES_PKHASH], &G2in[2][BYTES_PKHASH], &G2in[3][BYTES_PKHASH]
    };
    uint8_t *salt[4] = {
        &G2in[0][BYTES_PKHASH + BYTES_MU],
        &G2in[1][BYTES_PKHASH + BYTES_MU],
        &G2in[2][BYTES_PKHASH + BYTES_MU],
        &G2in[3][BYTES_PKHASH + BYTES_MU]
    };
    uint8_t G2out[4][BYTES_SEED_SE + CRYPTO_BYTES];                    // contains secret data
    uint8_t *seedSE[4] = { G2out[0], G2out[1], G2out[2], G2out[3] };   // contains secret data
    uint8_t *k[4] = {                                                  // contains secret data
        &G2out[0][BYTES_SEED_SE],
        &G2out[1][BYTES_SEED_SE],
        &G2out[2][BYTES_SEED_SE],
        &G2out[3][BYTES_SEED_SE]
    };
    uint8_t Fin[4][CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];             // contains secret data via Fin_k
    uint8_t *Fin_ct[4] = { Fin[0], Fin[1], Fin[2], Fin[3] };
    uint8_t *Fin_k[4] = {                                              // contains secret data
        &Fin[0][CRYPTO_CIPHERTEXTBYTES],
        &Fin[1][CRYPTO_CIPHERTEXTBYTES],
        &Fin[2][CRYPTO_CIPHERTEXTBYTES],
        &Fin[3][CRYPTO_CIPHERTEXTBYTES]
    };
    uint8_t shake_input_seedSE[4][1 + BYTES_SEED_SE];                  // contains secret data

#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_UNDEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif

#ifdef ALLOC_MMAP
    memset(B, 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
    for (size_t n = 0; n < 4; n++) {
        memset(V[n], 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        memset(Bp[n], 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        memset(Sp[n], 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        memset(Ep[n], 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        memset(Epp[n], 0, PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
    }
#endif

    frodo_unpack(B, PARAMS_N*PARAMS_NBAR, pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, PARAMS_LOGQ);

    // pkh <- G_1(pk)
    shakex2(pkh[0], pkh[1], BYTES_PKHASH, pk, pk, CRYPTO_PUBLICKEYBYTES);
    shakex2(pkh[2], pkh[3], BYTES_PKHASH, pk, pk, CRYPTO_PUBLICKEYBYTES);

    for (size_t n = 0; n < 4; n++) {
        // Generate random mu and salt
        if (randombytes(mu[n], BYTES_MU + BYTES_SALT) != 0)
            return 1;
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_UNDEFINED(mu[n], BYTES_MU + BYTES_SALT);
#endif
    }

    // Compute (seedSE || k) = G_2(pkh || mu || salt)
    shakex2(G2out[0], G2out[1], BYTES_SEED_SE + CRYPTO_BYTES, G2in[0], G2in[1], BYTES_PKHASH + BYTES_MU + BYTES_SALT);
    shakex2(G2out[2], G2out[3], BYTES_SEED_SE + CRYPTO_BYTES, G2in[2], G2in[3], BYTES_PKHASH + BYTES_MU + BYTES_SALT);

    for (size_t n = 0; n < 4; n++) {
        // Generate Sp and Ep, and compute Bp = Sp*A + Ep. Generate A on-the-fly
        shake_input_seedSE[n][0] = 0x96;
        memcpy(&shake_input_seedSE[n][1], seedSE[n], BYTES_SEED_SE);
    }

#ifdef ALLOC_MMAP
    keccakx2_state state[2];

    shakex2_absorb(&state[0], shake_input_seedSE[0], shake_input_seedSE[1], 1 + BYTES_SEED_SE);
    shakex2_absorb(&state[1], shake_input_seedSE[2], shake_input_seedSE[3], 1 + BYTES_SEED_SE);

    const size_t nblocks_Sp = (PARAMS_N*PARAMS_NBAR*sizeof(uint16_t) + SHAKE_RATE - 1)/SHAKE_RATE;
    const size_t offset_Ep = SHAKE_RATE*nblocks_Sp - PARAMS_N*PARAMS_NBAR*sizeof(uint16_t);
    const size_t nblocks_Ep = (PARAMS_N*PARAMS_NBAR*sizeof(uint16_t) - offset_Ep + SHAKE_RATE - 1)/SHAKE_RATE;
    const size_t offset_Epp = offset_Ep + SHAKE_RATE*nblocks_Ep - PARAMS_N*PARAMS_NBAR*sizeof(uint16_t);
    const size_t nblocks_Epp = (PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t) - offset_Epp + SHAKE_RATE - 1)/SHAKE_RATE;

    shakex2_squeezeblocks((uint8_t*)Sp[0], (uint8_t*)Sp[1], nblocks_Sp, &state[0]);
    shakex2_squeezeblocks((uint8_t*)Sp[2], (uint8_t*)Sp[3], nblocks_Sp, &state[1]);

    for (size_t n = 0; n < 4; n++) {
        memcpy(Ep[n], &Sp[n][PARAMS_N*PARAMS_NBAR], offset_Ep);
    }

    shakex2_squeezeblocks((uint8_t*)Ep[0] + offset_Ep, (uint8_t*)Ep[1] + offset_Ep, nblocks_Ep, &state[0]);
    shakex2_squeezeblocks((uint8_t*)Ep[2] + offset_Ep, (uint8_t*)Ep[3] + offset_Ep, nblocks_Ep, &state[1]);

    for (size_t n = 0; n < 4; n++) {
        memcpy(Epp[n], &Ep[n][PARAMS_N*PARAMS_NBAR], offset_Epp);
    }

    shakex2_squeezeblocks((uint8_t*)Epp[0] + offset_Epp, (uint8_t*)Epp[1] + offset_Epp, nblocks_Epp, &state[0]);
    shakex2_squeezeblocks((uint8_t*)Epp[2] + offset_Epp, (uint8_t*)Epp[3] + offset_Epp, nblocks_Epp, &state[1]);
#else
    shakex2(
        (uint8_t*)Sp[0], (uint8_t*)Sp[1], (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t),
        shake_input_seedSE[0], shake_input_seedSE[1], 1 + BYTES_SEED_SE
    );
    shakex2(
        (uint8_t*)Sp[2], (uint8_t*)Sp[3], (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t),
        shake_input_seedSE[2], shake_input_seedSE[3], 1 + BYTES_SEED_SE
    );
#endif

    for (size_t n = 0; n < 4; n++) {
        frodo_sample_n(Sp[n], PARAMS_N*PARAMS_NBAR);
        frodo_sample_n(Ep[n], PARAMS_N*PARAMS_NBAR);
    }

    frodo_mul_add_sa_plus_e_4x(Bp, (const uint16_t**)Sp, Ep, pk_seedA);
    
    for (size_t n = 0; n < 4; n++) {
        frodo_pack(ct_c1[n], (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8, Bp[n], PARAMS_N*PARAMS_NBAR, PARAMS_LOGQ);

        // Generate Epp, and compute V = Sp*B + Epp
        frodo_sample_n(Epp[n], PARAMS_NBAR*PARAMS_NBAR);
    }

    frodo_mul_add_sb_plus_e_4x(V, B, (const uint16_t**)Sp, (const uint16_t**)Epp);

    for (size_t n = 0; n < 4; n++) {
        // Encode mu, and compute C = V + enc(mu) (mod q)
        frodo_key_encode(C[n], (uint16_t*)mu[n]);
        frodo_add(C[n], V[n], C[n]);
        frodo_pack(ct_c2[n], (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8, C[n], PARAMS_NBAR*PARAMS_NBAR, PARAMS_LOGQ);

        // Append salt to ct and compute ss = F(ct_c1||ct_c2||salt||k)
        memcpy(&ct[n][CRYPTO_CIPHERTEXTBYTES - BYTES_SALT], salt[n], BYTES_SALT);
        memcpy(Fin_ct[n], ct[n], CRYPTO_CIPHERTEXTBYTES);
        memcpy(Fin_k[n], k[n], CRYPTO_BYTES);
    }

    shakex2(ss[0], ss[1], CRYPTO_BYTES, Fin[0], Fin[1], CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);
    shakex2(ss[2], ss[3], CRYPTO_BYTES, Fin[2], Fin[3], CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);    
    
    for (size_t n = 0; n < 4; n++) {
        // Cleanup:
        clear_bytes((uint8_t *)V[n], PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Sp[n], PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Ep[n], PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Epp[n], PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes(mu[n], BYTES_MU);
        clear_bytes(G2out[n], BYTES_SEED_SE + CRYPTO_BYTES);
        clear_bytes(Fin_k[n], CRYPTO_BYTES);
        clear_bytes(shake_input_seedSE[n], 1 + BYTES_SEED_SE);
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(mu, BYTES_MU);
#endif
    }

#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_DEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif
    return 0;
}

#ifdef ALLOC_MMAP
uint16_t *W_dec4x[4], *BBp_dec4x[4], *Sp_dec4x[4];

__attribute__((constructor)) static void alloc_dec4x(void) {
    for (int i = 0; i < 4; i++) {
        W_dec4x[i] = MEMORY_ALLOC(PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
        BBp_dec4x[i] = MEMORY_ALLOC(PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        Sp_dec4x[i] = MEMORY_ALLOC((2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR * sizeof(uint16_t));
    }
}
#endif

int crypto_kem_dec4x(unsigned char *ss0, const unsigned char *ct0, unsigned char *ss1, const unsigned char *ct1,
                     unsigned char *ss2, const unsigned char *ct2, unsigned char *ss3, const unsigned char *ct3,
                     const unsigned char *sk)
{ // FrodoKEM's key decapsulation
  // Inputs: ciphertext ct = ct_c1||ct_c2||salt                  (                              (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8 + BYTES_SALT bytes)
  //         secret key sk = sk_s||pk_seedA||pk_b||sk_S||sk_pkh  (CRYPTO_BYTES + BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + 2*PARAMS_N*PARAMS_NBAR + BYTES_PKHASH bytes)
  // Output: shared key ss                                       (CRYPTO_BYTES bytes)
    unsigned char* ss[4] = {ss0, ss1, ss2, ss3};
    const unsigned char* ct[4] = {ct0, ct1, ct2, ct3};
    uint16_t B[PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t Bp[4][PARAMS_N*PARAMS_NBAR] = {0};
    uint16_t C[4][PARAMS_NBAR*PARAMS_NBAR] = {0};
    uint16_t CC[4][PARAMS_NBAR*PARAMS_NBAR] = {0};
#ifdef ALLOC_MMAP
    uint16_t *W[4] = { W_dec4x[0], W_dec4x[1], W_dec4x[2], W_dec4x[3] }; // contains secret data
    uint16_t *BBp[4] = { BBp_dec4x[0], BBp_dec4x[1], BBp_dec4x[2], BBp_dec4x[3] };
    uint16_t *Sp[4] = { Sp_dec4x[0], Sp_dec4x[1], Sp_dec4x[2], Sp_dec4x[3] };
#else
    ALIGN_HEADER(32) uint16_t W0[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t W1[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t W2[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t W3[PARAMS_NBAR*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};

    uint16_t* W[4] = { W0, W1, W2, W3 };                                 // contains secret data

    ALIGN_HEADER(32) uint16_t BBp0[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t BBp1[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t BBp2[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) uint16_t BBp3[PARAMS_N*PARAMS_NBAR] ALIGN_FOOTER(32) = {0};

    uint16_t* BBp[4] = { BBp0, BBp1, BBp2, BBp3 };

    ALIGN_HEADER(128) uint16_t Sp0[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};
    ALIGN_HEADER(128) uint16_t Sp1[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};
    ALIGN_HEADER(128) uint16_t Sp2[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};
    ALIGN_HEADER(128) uint16_t Sp3[(2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR] ALIGN_FOOTER(128) = {0};

    uint16_t* Sp[4] = { Sp0, Sp1, Sp2, Sp3 };                          // contains secret data
#endif
    uint16_t *Ep[4] = {                                                // contains secret data
        (uint16_t *)&Sp[0][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[1][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[2][PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[3][PARAMS_N*PARAMS_NBAR]
    };
    uint16_t *Epp[4] = {                                               // contains secret data
        (uint16_t *)&Sp[0][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[1][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[2][2*PARAMS_N*PARAMS_NBAR],
        (uint16_t *)&Sp[3][2*PARAMS_N*PARAMS_NBAR]
    };
    const uint8_t *ct_c1[4] = { ct[0], ct[1], ct[2], ct[3] };
    const uint8_t *ct_c2[4] = {
        &ct[0][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[1][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[2][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8],
        &ct[3][(PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8]
    };
    const uint8_t *salt[4] = {
        &ct[0][CRYPTO_CIPHERTEXTBYTES - BYTES_SALT],
        &ct[1][CRYPTO_CIPHERTEXTBYTES - BYTES_SALT],
        &ct[2][CRYPTO_CIPHERTEXTBYTES - BYTES_SALT],
        &ct[3][CRYPTO_CIPHERTEXTBYTES - BYTES_SALT]
    };
    const uint8_t *sk_s = &sk[0];
    const uint8_t *sk_pk = &sk[CRYPTO_BYTES];
    const uint16_t *sk_S = (uint16_t *) &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES];
    uint16_t S[PARAMS_N * PARAMS_NBAR];                             // contains secret data
    const uint8_t *sk_pkh = &sk[CRYPTO_BYTES + CRYPTO_PUBLICKEYBYTES + 2*PARAMS_N*PARAMS_NBAR];
    const uint8_t *pk_seedA = &sk_pk[0];
    const uint8_t *pk_b = &sk_pk[BYTES_SEED_A];
    uint8_t G2in[4][BYTES_PKHASH + BYTES_MU + BYTES_SALT];             // contains secret data via muprime
    uint8_t *pkh[4] = { &G2in[0][0], &G2in[1][0], &G2in[2][0], &G2in[3][0] };
    uint8_t *muprime[4] = {                                            // contains secret data
        &G2in[0][BYTES_PKHASH], &G2in[1][BYTES_PKHASH], &G2in[2][BYTES_PKHASH], &G2in[3][BYTES_PKHASH]
    };
    uint8_t *G2in_salt[4] = {
        &G2in[0][BYTES_PKHASH + BYTES_MU],
        &G2in[1][BYTES_PKHASH + BYTES_MU],
        &G2in[2][BYTES_PKHASH + BYTES_MU],
        &G2in[3][BYTES_PKHASH + BYTES_MU]
    };
    uint8_t G2out[4][BYTES_SEED_SE + CRYPTO_BYTES];                    // contains secret data
    uint8_t *seedSEprime[4] = {                                        // contains secret data
        &G2out[0][0], &G2out[1][0], &G2out[2][0], &G2out[3][0]
    };
    uint8_t *kprime[4] = {                                             // contains secret data
        &G2out[0][BYTES_SEED_SE], &G2out[1][BYTES_SEED_SE], &G2out[2][BYTES_SEED_SE], &G2out[3][BYTES_SEED_SE]
    };
    uint8_t Fin[4][CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES];             // contains secret data via Fin_k
    uint8_t *Fin_ct[4] = { &Fin[0][0], &Fin[1][0], &Fin[2][0], &Fin[3][0] };
    uint8_t *Fin_k[4] = {                                              // contains secret data
        &Fin[0][CRYPTO_CIPHERTEXTBYTES],
        &Fin[1][CRYPTO_CIPHERTEXTBYTES],
        &Fin[2][CRYPTO_CIPHERTEXTBYTES],
        &Fin[3][CRYPTO_CIPHERTEXTBYTES]
    };
    uint8_t shake_input_seedSEprime[4][1 + BYTES_SEED_SE];              // contains secret data

#ifdef ALLOC_MMAP
    for (size_t n = 0; n < 4; n++) {
        memset(W[n], 0, PARAMS_NBAR*PARAMS_NBAR * sizeof(uint16_t));
        memset(BBp[n], 0, PARAMS_N*PARAMS_NBAR * sizeof(uint16_t));
        memset(Sp[n], 0, (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR * sizeof(uint16_t));
    }
#endif

#ifdef DO_VALGRIND_CHECK
    VALGRIND_MAKE_MEM_UNDEFINED(sk, CRYPTO_SECRETKEYBYTES);
    for (size_t n = 0; n < 4; n++) {
        VALGRIND_MAKE_MEM_UNDEFINED(ct[n], CRYPTO_CIPHERTEXTBYTES);
    }
#endif

    frodo_unpack(B, PARAMS_N*PARAMS_NBAR, pk_b, CRYPTO_PUBLICKEYBYTES - BYTES_SEED_A, PARAMS_LOGQ);

    for (size_t i = 0; i < PARAMS_N * PARAMS_NBAR; i++) {
        S[i] = LE_TO_UINT16(sk_S[i]);
    }

    for (size_t n = 0; n < 4; n++) {
        // Compute W = C - Bp*S (mod q), and decode the randomness mu
        frodo_unpack(Bp[n], PARAMS_N*PARAMS_NBAR, ct_c1[n], (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8, PARAMS_LOGQ);
        frodo_unpack(C[n], PARAMS_NBAR*PARAMS_NBAR, ct_c2[n], (PARAMS_LOGQ*PARAMS_NBAR*PARAMS_NBAR)/8, PARAMS_LOGQ);
        frodo_mul_bs(W[n], Bp[n], S);
        frodo_sub(W[n], C[n], W[n]);
        frodo_key_decode((uint16_t*)muprime[n], W[n]);

        // Generate (seedSE' || k') = G_2(pkh || mu' || salt)
        memcpy(pkh[n], sk_pkh, BYTES_PKHASH);
        memcpy(G2in_salt[n], salt[n], BYTES_SALT);
    }

    shakex2(G2out[0], G2out[1], BYTES_SEED_SE + CRYPTO_BYTES, G2in[0], G2in[1], BYTES_PKHASH + BYTES_MU + BYTES_SALT);
    shakex2(G2out[2], G2out[3], BYTES_SEED_SE + CRYPTO_BYTES, G2in[2], G2in[3], BYTES_PKHASH + BYTES_MU + BYTES_SALT);

    for (size_t n = 0; n < 4; n++) {
        // Generate Sp and Ep, and compute BBp = Sp*A + Ep. Generate A on-the-fly
        shake_input_seedSEprime[n][0] = 0x96;
        memcpy(&shake_input_seedSEprime[n][1], seedSEprime[n], BYTES_SEED_SE);
    }

    shakex2(
        (uint8_t*)Sp[0], (uint8_t*)Sp[1], (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t),
        shake_input_seedSEprime[0], shake_input_seedSEprime[1], 1 + BYTES_SEED_SE
    );

    shakex2(
        (uint8_t*)Sp[2], (uint8_t*)Sp[3], (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR*sizeof(uint16_t),
        shake_input_seedSEprime[2], shake_input_seedSEprime[3], 1 + BYTES_SEED_SE
    );

    for (size_t n = 0; n < 4; n++) {
        for (size_t i = 0; i < (2*PARAMS_N+PARAMS_NBAR)*PARAMS_NBAR; i++) {
            Sp[n][i] = LE_TO_UINT16(Sp[n][i]);
        }
        frodo_sample_n(Sp[n], PARAMS_N*PARAMS_NBAR);
        frodo_sample_n(Ep[n], PARAMS_N*PARAMS_NBAR);
    }

    frodo_mul_add_sa_plus_e_4x(BBp, (const uint16_t**)Sp, Ep, pk_seedA);

    for (size_t n = 0; n < 4; n++) {
        // Generate Epp, and compute W = Sp*B + Epp
        frodo_sample_n(Epp[n], PARAMS_NBAR*PARAMS_NBAR);
    }

    frodo_mul_add_sb_plus_e_4x(W, B, (const uint16_t**)Sp, (const uint16_t**)Epp);

    for (size_t n = 0; n < 4; n++) {
        // Encode mu, and compute CC = W + enc(mu') (mod q)
        frodo_key_encode(CC[n], (uint16_t*)muprime[n]);
        frodo_add(CC[n], W[n], CC[n]);

        // Prepare input to F
        memcpy(Fin_ct[n], ct[n], CRYPTO_CIPHERTEXTBYTES);

        // Reducing BBp modulo q
        for (int i = 0; i < PARAMS_N*PARAMS_NBAR; i++) BBp[n][i] = BBp[n][i] & ((1 << PARAMS_LOGQ)-1);

        // If (Bp == BBp & C == CC) then ss = F(ct || k'), else ss = F(ct || s)
        // Needs to avoid branching on secret data using constant-time implementation.
        int8_t selector = ct_verify(Bp[n], BBp[n], PARAMS_N*PARAMS_NBAR) | ct_verify(C[n], CC[n], PARAMS_NBAR*PARAMS_NBAR);
        // If (selector == 0) then load k' to do ss = F(ct || k'), else if (selector == -1) load s to do ss = F(ct || s)
        ct_select((uint8_t*)Fin_k[n], (uint8_t*)kprime[n], (uint8_t*)sk_s, CRYPTO_BYTES, selector);
    }

    shakex2(ss[0], ss[1], CRYPTO_BYTES, Fin[0], Fin[1], CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);
    shakex2(ss[2], ss[3], CRYPTO_BYTES, Fin[2], Fin[3], CRYPTO_CIPHERTEXTBYTES + CRYPTO_BYTES);
    
    for (size_t n = 0; n < 4; n++)
    {
        // Cleanup:
        clear_bytes((uint8_t *)W[n], PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Sp[n], PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Ep[n], PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes((uint8_t *)Epp[n], PARAMS_NBAR*PARAMS_NBAR*sizeof(uint16_t));
        clear_bytes(muprime[n], BYTES_MU);
        clear_bytes(G2out[n], BYTES_SEED_SE + CRYPTO_BYTES);
        clear_bytes(Fin_k[n], CRYPTO_BYTES);
        clear_bytes(shake_input_seedSEprime[n], 1 + BYTES_SEED_SE);
#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(ct[n], CRYPTO_CIPHERTEXTBYTES);
#endif
    }

#ifdef DO_VALGRIND_CHECK
        VALGRIND_MAKE_MEM_DEFINED(sk, CRYPTO_SECRETKEYBYTES);
#endif

    clear_bytes((uint8_t *)S, PARAMS_N*PARAMS_NBAR*sizeof(uint16_t));

    return 0;
}
