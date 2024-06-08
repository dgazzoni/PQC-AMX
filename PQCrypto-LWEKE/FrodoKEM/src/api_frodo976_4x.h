#ifndef _API_Frodo976_4x_H_
#define _API_Frodo976_4x_H_

#ifdef __cplusplus
extern "C" {
#endif

int crypto_kem_enc4x_Frodo976(unsigned char *ct0, unsigned char *ss0, unsigned char *ct1, unsigned char *ss1,
                              unsigned char *ct2, unsigned char *ss2, unsigned char *ct3, unsigned char *ss3,
                              const unsigned char *pk);

int crypto_kem_dec4x_Frodo976(unsigned char *ss0, const unsigned char *ct0, unsigned char *ss1,
                              const unsigned char *ct1, unsigned char *ss2, const unsigned char *ct2,
                              unsigned char *ss3, const unsigned char *ct3, const unsigned char *sk);

#ifdef __cplusplus
}
#endif

#endif
