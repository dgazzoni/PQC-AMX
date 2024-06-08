/********************************************************************************************
* FrodoKEM: Learning with Errors Key Encapsulation
*
* Abstract: header for internal functions
*********************************************************************************************/

#ifndef _FRODO_MACRIFY_4X_H_
#define _FRODO_MACRIFY_4X_H_

#include <stddef.h>
#include <stdint.h>
#include "config.h"


int frodo_mul_add_sa_plus_e_4x(uint16_t *b[4], const uint16_t *s[4], uint16_t *e[4], const uint8_t *seed_A);
void frodo_mul_add_sb_plus_e_4x(uint16_t *out[4], const uint16_t *b, const uint16_t *s[4], const uint16_t *e[4]);

#endif
