#include <cstring>
extern "C" {
#include "aes.h"
}

#include <cstdlib>

#include "gtest/gtest.h"
#include "test.h"

extern "C" void ARMv8_AES128_load_schedule(const uint8_t *key, uint8_t *schedule);
extern "C" void ARMv8_AES128_ECB_enc_sch(const uint8_t *plaintext, const size_t plaintext_len, const uint8_t *schedule,
                                         uint8_t *ciphertext);
extern "C" void AES128_ECB_genA_rows(const uint8_t *schedule, uint16_t *A, int n_rows, int first_row, int n,
                                     int n_cols);
extern "C" void AES128_ECB_genA_32cols(const uint8_t *schedule, uint16_t *A, int first_col, int n);

TEST(aes, FrodoKEM_impl_matches_FIPS197_C_1) {
    const uint8_t plaintext[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                                   0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};
    const uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                             0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t ciphertext_expected[16] = {0x69, 0xC4, 0xE0, 0xD8, 0x6A, 0x7B, 0x04, 0x30,
                                             0xD8, 0xCD, 0xB7, 0x80, 0x70, 0xB4, 0xC5, 0x5A};

    uint8_t schedule[16 * 11];
    uint8_t ciphertext[16];

    AES128_load_schedule(key, schedule);
    AES128_ECB_enc_sch(plaintext, 16, schedule, ciphertext);

    EXPECT_TRUE(ArraysMatch(ciphertext_expected, ciphertext));
}

TEST(aes, FrodoKEM_schedule_ARMv8_encrypt_matches_FIPS197_C_1_single_block) {
    const uint8_t plaintext[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                                   0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};
    const uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                             0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t ciphertext_expected[16] = {0x69, 0xC4, 0xE0, 0xD8, 0x6A, 0x7B, 0x04, 0x30,
                                             0xD8, 0xCD, 0xB7, 0x80, 0x70, 0xB4, 0xC5, 0x5A};

    uint8_t schedule[16 * 11];
    uint8_t ciphertext[16];

    AES128_load_schedule(key, schedule);
    AES128_ECB_enc_sch(plaintext, 16, schedule, ciphertext);

    EXPECT_TRUE(ArraysMatch(ciphertext_expected, ciphertext));
}

TEST(aes, ARMv8_impl_matches_FIPS197_C_1_single_block) {
    const uint8_t plaintext[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
                                   0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF};
    const uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                             0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    const uint8_t ciphertext_expected[16] = {0x69, 0xC4, 0xE0, 0xD8, 0x6A, 0x7B, 0x04, 0x30,
                                             0xD8, 0xCD, 0xB7, 0x80, 0x70, 0xB4, 0xC5, 0x5A};

    uint8_t schedule[16 * 11];
    uint8_t ciphertext[16];

    ARMv8_AES128_load_schedule(key, schedule);
    ARMv8_AES128_ECB_enc_sch(plaintext, 16, schedule, ciphertext);

    EXPECT_TRUE(ArraysMatch(ciphertext_expected, ciphertext));
}

static void test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(size_t n) {
    uint8_t key[16];
    uint8_t schedule[16 * 11], schedule_ARMv8[16 * 11];

    uint8_t *plaintext = new uint8_t[16 * n];
    uint8_t *ciphertext_expected = new uint8_t[16 * n];
    uint8_t *ciphertext = new uint8_t[16 * n];

    srand(0);

    for (size_t i = 0; i < 16 * n; i++) {
        plaintext[i] = rand() % 256;
    }

    for (size_t i = 0; i < 16; i++) {
        key[i] = rand() % 256;
    }

    AES128_load_schedule(key, schedule);
    AES128_ECB_enc_sch(plaintext, 16 * n, schedule, ciphertext_expected);

    ARMv8_AES128_load_schedule(key, schedule_ARMv8);
    ARMv8_AES128_ECB_enc_sch(plaintext, 16 * n, schedule_ARMv8, ciphertext);

    EXPECT_TRUE(ArraysMatch(ciphertext_expected, ciphertext, 16 * n));

    delete[] plaintext;
    delete[] ciphertext_expected;
    delete[] ciphertext;
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_two_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(2);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_three_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(3);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_four_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(4);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_five_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(5);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_six_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(6);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_seven_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(7);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_eight_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(8);
}

TEST(aes, ARMv8_impl_matches_FrodoKEM_impl_nine_blocks) {
    test_ARMv8_impl_matches_FrodoKEM_impl_n_blocks(9);
}

static void test_AES128_ECB_genA_rows(size_t rows, size_t cols, size_t n, size_t first_row) {
    uint8_t key[16];
    uint8_t schedule[16 * 11];

    uint16_t *Ainit = new uint16_t[rows * cols];
    uint16_t *A_expected = new uint16_t[rows * cols];
    uint16_t *A = new uint16_t[rows * cols];

    memset(Ainit, 0, rows * cols * sizeof(uint16_t));
    memset(A_expected, 0, rows * cols * sizeof(uint16_t));
    memset(A, 0, rows * cols * sizeof(uint16_t));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < n; j += 8) {
            Ainit[i * cols + j + 1] = j;
        }
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < n; j += 8) {
            Ainit[i * cols + j] = first_row + i;
        }
    }

    srand(1337);

    for (size_t i = 0; i < 16; i++) {
        key[i] = rand() % 256;
    }

    AES128_load_schedule(key, schedule);

    for (size_t i = 0; i < rows; i++) {
        AES128_ECB_enc_sch((uint8_t *)(Ainit + i * cols), n * sizeof(uint16_t), schedule,
                           (uint8_t *)(A_expected + i * cols));
    }

    AES128_ECB_genA_rows(schedule, A, rows, first_row, n, cols);

    EXPECT_TRUE(ArraysMatch(A_expected, A, rows * cols));

    delete[] Ainit;
    delete[] A_expected;
    delete[] A;
}

TEST(aes, AES128_ECB_genA_rows_4x640_firstrow_0) {
    test_AES128_ECB_genA_rows(4, 640, 640, 0);
}

TEST(aes, AES128_ECB_genA_rows_4x640_firstrow_636) {
    test_AES128_ECB_genA_rows(4, 640, 640, 636);
}

TEST(aes, AES128_ECB_genA_rows_8x640_first_row_0) {
    test_AES128_ECB_genA_rows(8, 640, 640, 0);
}

TEST(aes, AES128_ECB_genA_rows_8x640_first_row_632) {
    test_AES128_ECB_genA_rows(8, 640, 640, 632);
}

TEST(aes, AES128_ECB_genA_rows_32x640_first_row_0) {
    test_AES128_ECB_genA_rows(32, 640, 640, 0);
}

TEST(aes, AES128_ECB_genA_rows_32x640_first_row_608) {
    test_AES128_ECB_genA_rows(32, 640, 640, 608);
}

TEST(aes, AES128_ECB_genA_rows_4x976_first_row_0) {
    test_AES128_ECB_genA_rows(4, 976, 976, 0);
}

TEST(aes, AES128_ECB_genA_rows_4x976_first_row_972) {
    test_AES128_ECB_genA_rows(4, 976, 976, 972);
}

TEST(aes, AES128_ECB_genA_rows_8x976_first_row_0) {
    test_AES128_ECB_genA_rows(8, 976, 976, 0);
}

TEST(aes, AES128_ECB_genA_rows_8x976_first_row_968) {
    test_AES128_ECB_genA_rows(8, 976, 976, 968);
}

TEST(aes, AES128_ECB_genA_rows_4x976_992_first_row_0) {
    test_AES128_ECB_genA_rows(4, 992, 976, 0);
}

TEST(aes, AES128_ECB_genA_rows_4x976_992_first_row_972) {
    test_AES128_ECB_genA_rows(4, 992, 976, 972);
}

TEST(aes, AES128_ECB_genA_rows_8x976_992_first_row_0) {
    test_AES128_ECB_genA_rows(8, 992, 976, 0);
}

TEST(aes, AES128_ECB_genA_rows_8x976_992_first_row_968) {
    test_AES128_ECB_genA_rows(8, 992, 976, 968);
}

TEST(aes, AES128_ECB_genA_rows_32x976_992_first_row_0) {
    test_AES128_ECB_genA_rows(32, 992, 976, 0);
}

TEST(aes, AES128_ECB_genA_rows_16x976_992_first_row_960) {
    test_AES128_ECB_genA_rows(16, 992, 976, 960);
}

TEST(aes, AES128_ECB_genA_rows_4x1344_first_row_0) {
    test_AES128_ECB_genA_rows(4, 1344, 1344, 0);
}

TEST(aes, AES128_ECB_genA_rows_4x1344_first_row_1340) {
    test_AES128_ECB_genA_rows(4, 1344, 1344, 1340);
}

TEST(aes, AES128_ECB_genA_rows_8x1344_first_row_0) {
    test_AES128_ECB_genA_rows(8, 1344, 1344, 0);
}

TEST(aes, AES128_ECB_genA_rows_8x1344_first_row_1336) {
    test_AES128_ECB_genA_rows(8, 1344, 1344, 1336);
}

TEST(aes, AES128_ECB_genA_rows_32x1344_first_row_0) {
    test_AES128_ECB_genA_rows(32, 1344, 1344, 0);
}

TEST(aes, AES128_ECB_genA_rows_32x1344_first_row_1312) {
    test_AES128_ECB_genA_rows(32, 1344, 1344, 1312);
}

static void test_AES128_ECB_genA_32cols(size_t rows, int first_col) {
    uint8_t key[16];
    uint8_t schedule[16 * 11];

    uint16_t *Ainit = new uint16_t[rows * 32];
    uint16_t *A_expected = new uint16_t[rows * 32];
    uint16_t *A = new uint16_t[rows * 32];

    memset(Ainit, 0, rows * 32 * sizeof(uint16_t));
    memset(A_expected, 0, rows * 32 * sizeof(uint16_t));
    memset(A, 0, rows * 32 * sizeof(uint16_t));

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < 32; j += 8) {
            Ainit[i * 32 + j + 1] = first_col + j;
        }
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < 32; j += 8) {
            Ainit[i * 32 + j] = i;
        }
    }

    srand(1337);

    for (size_t i = 0; i < 16; i++) {
        key[i] = rand() % 256;
    }

    AES128_load_schedule(key, schedule);

    for (size_t i = 0; i < rows; i++) {
        AES128_ECB_enc_sch((uint8_t *)(Ainit + i * 32), 32 * sizeof(uint16_t), schedule,
                           (uint8_t *)(A_expected + i * 32));
    }

    AES128_ECB_genA_32cols(schedule, A, first_col, rows);

    EXPECT_TRUE(ArraysMatch(A_expected, A, rows * 32));

    delete[] Ainit;
    delete[] A_expected;
    delete[] A;
}

TEST(aes, AES128_ECB_genA_32cols_640x32_first_col_0) {
    test_AES128_ECB_genA_32cols(640, 0);
}

TEST(aes, AES128_ECB_genA_32cols_640x32_first_col_608) {
    test_AES128_ECB_genA_32cols(640, 608);
}

TEST(aes, AES128_ECB_genA_32cols_976x32_first_col_0) {
    test_AES128_ECB_genA_32cols(976, 0);
}

TEST(aes, AES128_ECB_genA_32cols_976x32_first_col_960) {
    test_AES128_ECB_genA_32cols(640, 960);
}

TEST(aes, AES128_ECB_genA_32cols_1344x32_first_col_0) {
    test_AES128_ECB_genA_32cols(1344, 0);
}

TEST(aes, AES128_ECB_genA_32cols_1344x32_first_col_1312) {
    test_AES128_ECB_genA_32cols(1344, 1312);
}
