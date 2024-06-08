extern "C" {
#include "rng.h"
}

#include "gtest/gtest.h"
#include "test.h"

TEST(TEST_NAME, test_kem4x_encaps4x1_decaps4x1) {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char ct[4][CRYPTO_CIPHERTEXTBYTES];
    unsigned char ss[4][CRYPTO_BYTES];
    unsigned char ss2[4][CRYPTO_BYTES];
    unsigned char entropy_input[48];

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    crypto_kem_keypair(pk, sk);
    crypto_kem_enc4x(ct[0], ss[0], ct[1], ss[1], ct[2], ss[2], ct[3], ss[3], pk);
    crypto_kem_dec4x(ss2[0], ct[0], ss2[1], ct[1], ss2[2], ct[2], ss2[3], ct[3], sk);

    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(ArraysMatch(ss[i], ss2[i]));
    }
}

TEST(TEST_NAME, test_kem4x_encaps1x4_decaps4x1) {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char ct[4][CRYPTO_CIPHERTEXTBYTES];
    unsigned char ss[4][CRYPTO_BYTES];
    unsigned char ss2[4][CRYPTO_BYTES];
    unsigned char entropy_input[48];

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    crypto_kem_keypair(pk, sk);

    for (int i = 0; i < 4; i++) {
        crypto_kem_enc(ct[i], ss[i], pk);
    }

    crypto_kem_dec4x(ss2[0], ct[0], ss2[1], ct[1], ss2[2], ct[2], ss2[3], ct[3], sk);

    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(ArraysMatch(ss[i], ss2[i]));
    }
}

TEST(TEST_NAME, test_kem4x_encaps4x1_decaps1x4) {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char ct[4][CRYPTO_CIPHERTEXTBYTES];
    unsigned char ss[4][CRYPTO_BYTES];
    unsigned char ss2[4][CRYPTO_BYTES];
    unsigned char entropy_input[48];

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    crypto_kem_keypair(pk, sk);

    crypto_kem_enc4x(ct[0], ss[0], ct[1], ss[1], ct[2], ss[2], ct[3], ss[3], pk);

    for (int i = 0; i < 4; i++) {
        crypto_kem_dec(ss2[i], ct[i], sk);
    }

    for (int i = 0; i < 4; i++) {
        EXPECT_TRUE(ArraysMatch(ss[i], ss2[i]));
    }
}
