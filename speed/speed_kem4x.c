
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "rng.h"

#ifndef NTESTS
#error Please #define NTESTS to the desired number of repetitions for the benchmark
#endif

uint64_t time0, time1;
uint64_t cycles[NTESTS];

#ifdef __APPLE__

#include "m1cycles.h"
#define SETUP_COUNTER() \
    {                   \
        (void)cycles;   \
        setup_rdtsc();  \
    }
#define CYCLE_TYPE "%lld"
#define GET_TIME rdtsc()

#else

#include "hal.h"
#define SETUP_COUNTER() \
    {}
#define CYCLE_TYPE "%ld"
#define GET_TIME hal_get_time()

#endif

#undef __AVERAGE__
#define __MEDIAN__

#ifdef __AVERAGE__

#define LOOP_INIT(__clock0, __clock1) \
    { __clock0 = GET_TIME; }
#define LOOP_TAIL(__f_string, records, __clock0, __clock1)  \
    {                                                       \
        __clock1 = GET_TIME;                                \
        printf(__f_string, (__clock1 - __clock0) / NTESTS); \
    }
#define BODY_INIT(__clock0, __clock1) \
    {}
#define BODY_TAIL(records, __clock0, __clock1) \
    {}

#elif defined(__MEDIAN__)

static int cmp_uint64(const void *a, const void *b) {
    return ((*((const uint64_t *)a)) - ((*((const uint64_t *)b))));
}

#define LOOP_INIT(__clock0, __clock1) \
    {}
#define LOOP_TAIL(__f_string, records, __clock0, __clock1)    \
    {                                                         \
        qsort(records, sizeof(uint64_t), NTESTS, cmp_uint64); \
        printf(__f_string, records[NTESTS >> 1]);             \
    }
#define BODY_INIT(__clock0, __clock1) \
    { __clock0 = GET_TIME; }
#define BODY_TAIL(records, __clock0, __clock1) \
    {                                          \
        __clock1 = GET_TIME;                   \
        records[i] = __clock1 - __clock0;      \
    }

#endif

#define WRAP_FUNC(__f_string, records, __clock0, __clock1, func) \
    {                                                            \
        LOOP_INIT(__clock0, __clock1);                           \
        for (size_t i = 0; i < NTESTS; i++) {                    \
            BODY_INIT(__clock0, __clock1);                       \
            func;                                                \
            BODY_TAIL(records, __clock0, __clock1);              \
        }                                                        \
        LOOP_TAIL(__f_string, records, __clock0, __clock1);      \
    }

int main() {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES] = {0};
    unsigned char sk[CRYPTO_SECRETKEYBYTES] = {0};
    unsigned char ct[4][CRYPTO_CIPHERTEXTBYTES] = {0};
    unsigned char key_a[4][CRYPTO_BYTES] = {0};
    unsigned char key_b[4][CRYPTO_BYTES] = {0};
    unsigned char entropy_input[48] = {0};

    for (int i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    randombytes_init(entropy_input, NULL, 256);

    crypto_kem_keypair(pk, sk);

    SETUP_COUNTER();

    WRAP_FUNC("crypto_kem_enc4x: " CYCLE_TYPE "\n", cycles, time0, time1,
              crypto_kem_enc4x(ct[0], key_b[0], ct[1], key_b[1], ct[2], key_b[2], ct[3], key_b[3], pk));
    WRAP_FUNC("crypto_kem_dec4x: " CYCLE_TYPE "\n", cycles, time0, time1,
              crypto_kem_dec4x(key_a[0], ct[0], key_a[1], ct[1], key_a[2], ct[2], key_a[3], ct[3], sk));

    return 0;
}