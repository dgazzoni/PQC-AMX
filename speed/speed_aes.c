
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "aes.h"
#include "speed.h"

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

#define WRAP_FUNC(__f_string, records, __clock0, __clock1, func, dno_target) \
    {                                                                        \
        LOOP_INIT(__clock0, __clock1);                                       \
        for (size_t i = 0; i < NTESTS; i++) {                                \
            BODY_INIT(__clock0, __clock1);                                   \
            func;                                                            \
            DoNotOptimize(dno_target);                                       \
            BODY_TAIL(records, __clock0, __clock1);                          \
        }                                                                    \
        LOOP_TAIL(__f_string, records, __clock0, __clock1);                  \
    }

int main() {
    uint8_t plaintext[65536], ciphertext[65536], key[16] = {0}, schedule[16 * 11];

    AES128_load_schedule(key, schedule);

    SETUP_COUNTER();

    WRAP_FUNC("AES128_ECB_enc_sch: " CYCLE_TYPE "\n", cycles, time0, time1,
              AES128_ECB_enc_sch(plaintext, sizeof(plaintext), schedule, ciphertext), ciphertext);

    return 0;
}
