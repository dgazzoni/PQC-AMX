#include "api_frodo640.h"
#include "api_frodo640_4x.h"
#include "frodo_macrify.h"
#include "frodo_macrify_4x.h"

// Parameters for "FrodoKEM-640"
#define PARAMS_N 640
#define PARAMS_NBAR 8
#define PARAMS_LOGQ 15
#define PARAMS_Q (1 << PARAMS_LOGQ)
#define PARAMS_EXTRACTED_BITS 2
#define PARAMS_STRIPE_STEP 8
#define PARAMS_PARALLEL 4
#define BYTES_SEED_A 16
#define BYTES_MU (PARAMS_EXTRACTED_BITS*PARAMS_NBAR*PARAMS_NBAR)/8
#define BYTES_SALT (2*CRYPTO_BYTES)
#define BYTES_SEED_SE (2*CRYPTO_BYTES)
#define BYTES_PKHASH CRYPTO_BYTES

#define shakex2     shake128x2
#define SHAKE_RATE  SHAKE128_RATE

#define crypto_kem_enc4x              crypto_kem_enc4x_Frodo640
#define crypto_kem_dec4x              crypto_kem_dec4x_Frodo640

#include "kem4x.c"
#if defined(AMX)
#if defined(USE_REFERENCE)
#include "frodo_macrify_reference_amx_4x.c"
#else
#include "frodo_macrify_amx_4x.c"
#endif
#else
#if defined(USE_REFERENCE)
#include "frodo_macrify_reference_4x.c"
#else
#include "frodo_macrify_4x.c"
#endif
#endif
