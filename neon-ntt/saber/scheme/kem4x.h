/*=============================================================================
This file has been adapted from the implementation
(available at, Public Domain https://github.com/KULeuven-COSIC/SABER)
of "Saber: Module-LWR based key exchange, CPA-secure encryption and CCA-secure KEM"
by : Jan-Pieter D'Anvers, Angshuman Karmakar, Sujoy Sinha Roy, and Frederik Vercauteren
Jose Maria Bermudo Mera, Michiel Van Beirendonck, Andrea Basso.
=============================================================================*/

#ifndef KEM4X_H
#define KEM4X_H

#include <stdint.h>

#include "kem.h"

#ifdef __cplusplus
extern "C" {
#endif

int crypto_kem_keypair4x(unsigned char *pk0, unsigned char *sk0, unsigned char *pk1, unsigned char *sk1,
                         unsigned char *pk2, unsigned char *sk2, unsigned char *pk3, unsigned char *sk3);

int crypto_kem_enc4x(unsigned char *c0, unsigned char *k0, const unsigned char *pk0, unsigned char *c1,
                     unsigned char *k1, const unsigned char *pk1, unsigned char *c2, unsigned char *k2,
                     const unsigned char *pk2, unsigned char *c3, unsigned char *k3, const unsigned char *pk3);

int crypto_kem_dec4x(unsigned char *k0, const unsigned char *c0, const unsigned char *sk0, unsigned char *k1,
                     const unsigned char *c1, const unsigned char *sk1, unsigned char *k2, const unsigned char *c2,
                     const unsigned char *sk2, unsigned char *k3, const unsigned char *c3, const unsigned char *sk3);

#ifdef __cplusplus
}
#endif

#endif
