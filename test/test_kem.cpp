extern "C" {
#include "rng.h"
}

#include "gtest/gtest.h"
#include "test.h"

TEST(TEST_NAME, test_kem) {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
    unsigned char ss[CRYPTO_BYTES];
    unsigned char ss2[CRYPTO_BYTES];
    unsigned char entropy_input[48];

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    crypto_kem_keypair(pk, sk);
    crypto_kem_enc(ct, ss, pk);
    crypto_kem_dec(ss2, ct, sk);

    EXPECT_TRUE(ArraysMatch(ss, ss2));
}
