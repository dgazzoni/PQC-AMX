#include "kem.h"

int crypto_kem_keypair4x(unsigned char *pk0, unsigned char *sk0, unsigned char *pk1, unsigned char *sk1,
                         unsigned char *pk2, unsigned char *sk2, unsigned char *pk3, unsigned char *sk3) {
    crypto_kem_keypair(pk0, sk0);
    crypto_kem_keypair(pk1, sk1);
    crypto_kem_keypair(pk2, sk2);
    crypto_kem_keypair(pk3, sk3);

    return 0;
}

int crypto_kem_enc4x(unsigned char *c0, unsigned char *k0, const unsigned char *pk0, unsigned char *c1,
                     unsigned char *k1, const unsigned char *pk1, unsigned char *c2, unsigned char *k2,
                     const unsigned char *pk2, unsigned char *c3, unsigned char *k3, const unsigned char *pk3) {
    crypto_kem_enc(c0, k0, pk0);
    crypto_kem_enc(c1, k1, pk1);
    crypto_kem_enc(c2, k2, pk2);
    crypto_kem_enc(c3, k3, pk3);

    return 0;
}

int crypto_kem_dec4x(unsigned char *k0, const unsigned char *c0, const unsigned char *sk0, unsigned char *k1,
                     const unsigned char *c1, const unsigned char *sk1, unsigned char *k2, const unsigned char *c2,
                     const unsigned char *sk2, unsigned char *k3, const unsigned char *c3, const unsigned char *sk3) {
    crypto_kem_dec(k0, c0, sk0);
    crypto_kem_dec(k1, c1, sk1);
    crypto_kem_dec(k2, c2, sk2);
    crypto_kem_dec(k3, c3, sk3);

    return 0;
}
