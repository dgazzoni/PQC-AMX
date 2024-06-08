#include <arm_neon.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>

#define ROTR32(x, n) ((x << (32 - n)) | (x >> n))

typedef union {
    uint8_t u8[11][16];
    uint32_t u32[11][4];
} subkeys_t;

static inline uint32_t AES_sbox_x4(uint32_t in) {
    uint8x16_t sbox_val = vreinterpretq_u8_u32(vdupq_n_u32(in));
    sbox_val = vaeseq_u8(sbox_val, vdupq_n_u8(0));

    return vgetq_lane_u32(vreinterpretq_u32_u8(sbox_val), 0);
}

void AES128_load_schedule(const uint8_t *key, uint8_t *subkeys) {
    subkeys_t *sk = (subkeys_t *)subkeys;
    uint8_t rcon = 1;
    uint32_t s;
    int i, j;

    memcpy(subkeys, key, 16 * sizeof(uint8_t));

    for (i = 1; i <= 10; i++) {
        s = AES_sbox_x4(sk->u32[i - 1][3]);
        sk->u32[i][0] = ROTR32(s, 8) ^ rcon ^ sk->u32[i - 1][0];

        for (j = 1; j < 4; j++) {
            sk->u32[i][j] = sk->u32[i][j - 1] ^ sk->u32[i - 1][j];
        }

        rcon = (rcon << 1) ^ ((rcon >> 7) * 0x11b);
    }
}

#define AES128_ECB_XWAYS(ways, vsubkeys, ctr, out)           \
    do {                                                     \
        uint8x16_t state[ways];                              \
                                                             \
        for (int j = 0; j < ways; j++) {                     \
            state[j] = vaeseq_u8(ctr[j], vsubkeys[0]);       \
            state[j] = vaesmcq_u8(state[j]);                 \
        }                                                    \
                                                             \
        for (int i = 1; i < 9; i++) {                        \
            for (int j = 0; j < ways; j++) {                 \
                state[j] = vaeseq_u8(state[j], vsubkeys[i]); \
                state[j] = vaesmcq_u8(state[j]);             \
            }                                                \
        }                                                    \
                                                             \
        for (int j = 0; j < ways; j++) {                     \
            state[j] = vaeseq_u8(state[j], vsubkeys[9]);     \
            state[j] = veorq_u8(state[j], vsubkeys[10]);     \
            vst1q_u8(out + j * 16, state[j]);                \
        }                                                    \
    }                                                        \
    while (0);

static void AES128_ECB(uint8x16_t vsubkeys[11], uint8x16_t in, unsigned char *out) {
    AES128_ECB_XWAYS(1, vsubkeys, (&in), out);
}

static void AES128_ECB_x8(uint8x16_t vsubkeys[11], uint8x16_t in[8], unsigned char *out) {
    AES128_ECB_XWAYS(8, vsubkeys, in, out);
}

// TODO: conjecturally, the optimal number of ways is 12, but this requires very precise instruction scheduling and
// alignment, which is better done in assembly language; try this in the future
void AES128_ECB_enc_sch(const uint8_t *plaintext, const size_t plaintext_len, const uint8_t *schedule,
                        uint8_t *ciphertext) {
    size_t block;
    uint8x16_t vsubkeys[11], vplaintext[8];

    assert(plaintext_len % 16 == 0);

    for (int i = 0; i < 11; i++) {
        vsubkeys[i] = vld1q_u8(schedule + i * 16);
    }

    for (block = 0; block < plaintext_len / 128; block++) {
        for (int i = 0; i < 8; i++) {
            vplaintext[i] = vld1q_u8(&plaintext[128 * block + 16 * i]);
        }
        AES128_ECB_x8(vsubkeys, vplaintext, &ciphertext[128 * block]);
    }

    for (block = 8 * block; block < plaintext_len / 16; block++) {
        vplaintext[0] = vld1q_u8(&plaintext[16 * block]);
        AES128_ECB(vsubkeys, vplaintext[0], &ciphertext[16 * block]);
    }
}

void AES128_ECB_genA_rows(const uint8_t *schedule, uint16_t *A, int n_rows, int first_row, int n, int n_cols) {
    int i, j, k;
    uint32x4_t vi = vsetq_lane_u32(first_row, vdupq_n_u32(0), 0);
    uint32x4_t v1i = vsetq_lane_u32(1, vdupq_n_u32(0), 0);
    uint32x4_t v8j = vsetq_lane_u32(0x80000, vdupq_n_u32(0), 0);
    uint32x4_t vAinit = vi;
    uint8x16_t vsubkeys[11], vplaintext[8];

    for (i = 0; i < 11; i++) {
        vsubkeys[i] = vld1q_u8(schedule + i * 16);
    }

    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < 64 * (n / 64); j += 64) {
            for (k = 0; k < 8; k++) {
                vplaintext[k] = vreinterpretq_u8_u32(vAinit);
                vAinit = vaddq_u32(vAinit, v8j);
            }
            AES128_ECB_x8(vsubkeys, vplaintext, (uint8_t *)&A[i * n_cols + j]);
        }

        for (; j < n; j += 8) {
            vplaintext[0] = vreinterpretq_u8_u32(vAinit);
            vAinit = vaddq_u32(vAinit, v8j);
            AES128_ECB(vsubkeys, vplaintext[0], (uint8_t *)&A[i * n_cols + j]);
        }

        vi = vaddq_u32(vi, v1i);
        vAinit = vi;
    }
}

void AES128_ECB_genA_32cols(const uint8_t *schedule, uint16_t *A, int first_col, int n) {
    // Valid values for cols: 16 or 32
    int i, j;
    uint32x4_t vi = vsetq_lane_u32(first_col << 16, vdupq_n_u32(0), 0);
    uint32x4_t v1i = vsetq_lane_u32(1, vdupq_n_u32(0), 0);
    uint32x4_t v8j = vsetq_lane_u32(0x80000, vdupq_n_u32(0), 0);
    uint32x4_t vAinit = vi;
    uint8x16_t vsubkeys[11], vplaintext[8];

    for (i = 0; i < 11; i++) {
        vsubkeys[i] = vld1q_u8(schedule + i * 16);
    }

    for (i = 0; i < n; i += 2) {
        for (j = 0; j < 4; j++) {
            vplaintext[j] = vreinterpretq_u8_u32(vAinit);
            vAinit = vaddq_u32(vAinit, v8j);
        }

        vi = vaddq_u32(vi, v1i);
        vAinit = vi;

        for (j = 4; j < 8; j++) {
            vplaintext[j] = vreinterpretq_u8_u32(vAinit);
            vAinit = vaddq_u32(vAinit, v8j);
        }

        vi = vaddq_u32(vi, v1i);
        vAinit = vi;

        AES128_ECB_x8(vsubkeys, vplaintext, (uint8_t *)&A[i * 32]);
    }
}

void AES128_free_schedule(uint8_t *schedule) {
    memset(schedule, 0, 16 * 11);
}
