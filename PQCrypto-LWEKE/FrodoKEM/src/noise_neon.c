/********************************************************************************************
 * FrodoKEM: Learning with Errors Key Encapsulation
 *
 * Abstract: noise sampling functions
 *********************************************************************************************/

#include <arm_neon.h>

#include "../../common/sha3/fips202.h"

void frodo_sample_n(uint16_t *s, const size_t n) {
    // Fills vector s with n samples from the noise distribution which requires 16 bits to sample.
    // The distribution is specified by its CDF.
    // Input: pseudo-random values (2*n bytes) passed in s. The input is overwritten by the output.
    unsigned int i, j;
    uint16x8_t vcdf_table[CDF_TABLE_LEN - 1];

    for (i = 0; i < (unsigned int)(CDF_TABLE_LEN - 1); i++) {
        vcdf_table[i] = vdupq_n_u16(CDF_TABLE[i]);
    }

    for (i = 0; i < n; i += 8) {
        uint16x8_t vsample = vdupq_n_u16(0), vs;
        vs = vld1q_u16(&s[i]);
        uint16x8_t vprnd = vshrq_n_u16(vs, 1);  // Drop the least significant bit
        uint16x8_t vsign;
        asm volatile("cmtst.8h %0, %1, %2"
                     : "=w"(vsign)
                     : "w"(vs), "w"(vdupq_n_u16(0x1)));  // Pick the least significant bit

        // No need to compare with the last value.
        for (j = 0; j < (unsigned int)(CDF_TABLE_LEN - 1); j++) {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit
            // in 15 bits.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
            vsample = vsubq_u16(vsample, vcltq_u16(vcdf_table[j], vprnd));
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1
        asm volatile("bsl.16b %0, %1, %2" : "+w"(vsign) : "w"(vnegq_s16(vreinterpretq_s16_u16(vsample))), "w"(vsample));
        vst1q_u16(&s[i], vsign);
    }
}
