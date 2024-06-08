#ifndef AMX_H
#define AMX_H

#include "aarch64.h"

#define AMX_PTR(ptr) ((uint64_t)(ptr))

// loads
// Default arguments for load instructions: LDX_LOAD_SINGLE, LDY_LOAD_SINGLE, LDZ_LOAD_SINGLE

// ldx

#define LDX_REG(reg) ((uint64_t)(reg) << 56)

#define LDX_LOAD_SINGLE (0ULL << 62)
#define LDX_LOAD_PAIR (1ULL << 62)

// ldy

#define LDY_REG(reg) ((uint64_t)(reg) << 56)

#define LDY_LOAD_SINGLE (0ULL << 62)
#define LDY_LOAD_PAIR (1ULL << 62)

#ifdef CPU_M1
// M1 requires 128-byte aligned address for a load pair (https://github.com/corsix/amx/blob/main/ldst.md), which is
// being used to emulate a load quad. In general that's not the case, so we use load single.
#define AMX_LDX_QUAD(PTR, REG)                                              \
    do {                                                                    \
        for (int j = 0; j < 4; j++) {                                       \
            AMX_LDX(AMX_PTR(((uint8_t *)PTR + 64 * j)) | LDX_REG(REG + j)); \
        }                                                                   \
    }                                                                       \
    while (0)

#define AMX_LDY_QUAD(PTR, REG)                                              \
    do {                                                                    \
        for (int j = 0; j < 4; j++) {                                       \
            AMX_LDY(AMX_PTR(((uint8_t *)PTR + 64 * j)) | LDY_REG(REG + j)); \
        }                                                                   \
    }                                                                       \
    while (0)
#else
#define LDX_LOAD_QUAD (LDX_LOAD_PAIR | (1ULL << 60))
#define LDY_LOAD_QUAD (LDY_LOAD_PAIR | (1ULL << 60))

#define AMX_LDX_QUAD(PTR, REG)                                \
    do {                                                      \
        AMX_LDX(AMX_PTR(PTR) | LDX_REG(REG) | LDX_LOAD_QUAD); \
    }                                                         \
    while (0)

#define AMX_LDY_QUAD(PTR, REG)                                \
    do {                                                      \
        AMX_LDY(AMX_PTR(PTR) | LDY_REG(REG) | LDY_LOAD_QUAD); \
    }                                                         \
    while (0)
#endif

// ldz

#define LDZ_Z_ROW(row) ((uint64_t)(row) << 56)

#define LDZ_LOAD_SINGLE (0ULL << 62)
#define LDZ_LOAD_PAIR (1ULL << 62)

// stx

#define STX_REG(reg) ((uint64_t)(reg) << 56)

#define STX_STORE_SINGLE (0ULL << 62)
#define STX_STORE_PAIR (1ULL << 62)

// sty

#define STY_REG(reg) ((uint64_t)(reg) << 56)

#define STY_STORE_SINGLE (0ULL << 62)
#define STY_STORE_PAIR (1ULL << 62)

// stz
// Defalt arguments: STZ_STORE_SINGLE

#define STZ_Z_ROW(row) ((uint64_t)(row) << 56)

#define STZ_STORE_SINGLE (0ULL << 62)
#define STZ_STORE_PAIR (1ULL << 62)

// mac16
// Default arguments: MAC16_MATRIX, MAC16_X_NO_SKIP, MAC16_Y_NO_SKIP, MAC16_Z_NO_SKIP, MAC16_X_ENABLE_MODE(0),
//                    MAC16_X_ENABLE_VALUE(0), MAC16_MATRIX_Y_ENABLE_MODE(0), MAC16_MATRIX_Y_ENABLE_VALUE(0)

#define MAC16_MATRIX (0ULL << 63)
#define MAC16_VECTOR (1ULL << 63)

// mac16 common

#define MAC16_X_OFFSET(bytes) ((uint64_t)(bytes) << 10)
#define MAC16_X_REG(reg) (MAC16_X_OFFSET((reg)*64))

#define MAC16_Y_OFFSET(bytes) ((uint64_t)(bytes) << 0)
#define MAC16_Y_REG(reg) (MAC16_Y_OFFSET((reg)*64))

#define MAC16_Z_ROW(row) ((uint64_t)(row) << 20)

#define MAC16_Z_NO_SKIP (0ULL << 27)
#define MAC16_Z_SKIP (1ULL << 27)

#define MAC16_Y_NO_SKIP (0ULL << 28)
#define MAC16_Y_SKIP (1ULL << 28)

#define MAC16_X_NO_SKIP (0ULL << 29)
#define MAC16_X_SKIP (1ULL << 29)

#define MAC16_X_ENABLE_VALUE(value) ((uint64_t)(value) << 41)

#define MAC16_X_ENABLE_MODE(mode) ((uint64_t)(mode) << 46)

// mac16 matrix mode (bit 63 = 0)

#define MAC16_MATRIX_Y_ENABLE_VALUE(value) ((uint64_t)(value) << 32)

#define MAC16_MATRIX_Y_ENABLE_MODE(mode) ((uint64_t)(mode) << 37)

// vecint
// Default arguments: VECINT_ALU_MODE_Z_ADD_X_MUL_Y

#define VECINT_X_OFFSET(bytes) ((uint64_t)(bytes) << 10)
#define VECINT_X_REG(reg) (VECINT_X_OFFSET((reg)*64))

#define VECINT_Y_OFFSET(bytes) ((uint64_t)(bytes) << 0)
#define VECINT_Y_REG(reg) (VECINT_Y_OFFSET((reg)*64))

#define VECINT_Z_ROW(row) ((uint64_t)(row) << 20)

#define VECINT_ENABLE_VALUE(value) ((uint64_t)(value) << 32)

#define VECINT_ENABLE_MODE(mode) ((uint64_t)(mode) << 38)

#define VECINT_BROADCAST_MODE(mode) ((uint64_t)(mode) << 32)

#define VECINT_MULTIPLE_1 (0ULL << 31)
#define VECINT_MULTIPLE_2 ((1ULL << 31) | (0ULL << 25))
#define VECINT_MULTIPLE_4 ((1ULL << 31) | (1ULL << 25))

#define VECINT_ALU_MODE_Z_ADD_X_MUL_Y (0ULL << 47)
#define VECINT_ALU_MODE_Z_ADD_X_ADD_Y (2ULL << 47)
#define VECINT_ALU_MODE_Z_SUB_X_ADD_Y (3ULL << 47)
#define VECINT_ALU_MODE_Z_SHIFT_RIGHT_S (4ULL << 47)
#define VECINT_ALU_MODE_X_MUL_Y (10ULL << 47)

#define VECINT_RIGHT_SHIFT_AMOUNT(amount) ((uint64_t)(amount) << 58)

#define VECINT_Y_UNSIGNED (0ULL << 26)
#define VECINT_Y_SIGNED (1ULL << 26)

#define VECINT_X_UNSIGNED (0ULL << 63)
#define VECINT_X_SIGNED (1ULL << 63)

// extrh

// extrh not defined in aarch64.h
#define AMX_EXTRH(gpr) AMX_OP_GPR(8, gpr)

#define EXTRH_COPY_ONLY (0ULL << 26)
#define EXTRH_GENERIC (1ULL << 26)

// extrh copy-only mode (bit 26 = 0)

#define EXTRH_COPY_ONLY_OFFSET(offset) ((uint64_t)(offset) << 10)
#define EXTRH_COPY_ONLY_REG(reg) EXTRH_COPY_ONLY_OFFSET((reg)*64)

#define EXTRH_COPY_ONLY_Z_ROW(row) ((uint64_t)(row) << 20)

#define EXTRH_COPY_ONLY_LANE_WIDTH_16_BIT (2ULL << 28)

#define EXTRH_COPY_ONLY_ENABLE_VALUE(value) ((uint64_t)(value) << 41)

#define EXTRH_COPY_ONLY_ENABLE_MODE(mode) ((uint64_t)(mode) << 46)

// extrh generic mode (bit 26 = 1)
// Default arguments for extrh generic: EXTRH_DESTINATION_X

#define EXTRH_OFFSET(offset) ((uint64_t)(offset) << 0)
#define EXTRH_REG(reg) EXTRH_OFFSET((reg)*64)

#define EXTRH_DESTINATION_X (0ULL << 10)
#define EXTRH_DESTINATION_Y (1ULL << 10)

#define EXTRH_Z_ROW(row) ((uint64_t)(row) << 20)

#define EXTRH_MULTIPLE_1 (0ULL << 31)
#define EXTRH_MULTIPLE_2 ((1ULL << 31) | (0ULL << 25))
#define EXTRH_MULTIPLE_4 ((1ULL << 31) | (1ULL << 25))

#define EXTRH_LANE_WIDTH_16_BIT (1ULL << 11)

// extrv

// extrv not defined in aarch64.h
#define AMX_EXTRV(gpr) AMX_OP_GPR(9, gpr)

#define EXTRV_COPY_ONLY (0ULL << 26)
#define EXTRV_GENERIC (1ULL << 26)

// extrv copy-only mode (bit 26 = 0)

#define EXTRV_COPY_ONLY_OFFSET(offset) ((uint64_t)(offset) << 0)
#define EXTRV_COPY_ONLY_REG(reg) EXTRH_COPY_ONLY_OFFSET((reg)*64)

#define EXTRV_COPY_ONLY_Z_COLUMN(column) ((uint64_t)(column) << 20)

#define EXTRV_COPY_ONLY_LANE_WIDTH_16_BIT (2ULL << 28)

#define EXTRV_COPY_ONLY_ENABLE_VALUE(value) ((uint64_t)(value) << 32)

#define EXTRV_COPY_ONLY_ENABLE_MODE(mode) ((uint64_t)(mode) << 37)

// extrv generic mode (bit 26 = 1)
// Default arguments for extrh generic: EXTRV_DESTINATION_X

#define EXTRV_OFFSET(offset) ((uint64_t)(offset) << 0)
#define EXTRV_REG(reg) EXTRH_OFFSET((reg)*64)

#define EXTRV_DESTINATION_X (0ULL << 10)
#define EXTRV_DESTINATION_Y (1ULL << 10)

#define EXTRV_Z_COLUMN(column) ((uint64_t)(column) << 20)

#define EXTRV_ENABLE_VALUE(value) ((uint64_t)(value) << 32)

#define EXTRV_ENABLE_MODE(mode) ((uint64_t)(mode) << 38)

#define EXTRV_MULTIPLE_1 (0ULL << 31)
#define EXTRV_MULTIPLE_2 ((1ULL << 31) | (0ULL << 25))
#define EXTRV_MULTIPLE_4 ((1ULL << 31) | (1ULL << 25))

#define EXTRV_LANE_WIDTH_16_BIT (1ULL << 11)

// matint
#define MATINT_ALU_MODE_Z_ADD_X_MUL_Y (0ULL << 47)
#define MATINT_ALU_MODE_Z_SUB_X_MUL_Y (1ULL << 47)

#define MATINT_X_OFFSET(bytes) ((uint64_t)(bytes) << 10)
#define MATINT_X_REG(reg) (MATINT_X_OFFSET((reg)*64))

#define MATINT_Y_OFFSET(bytes) ((uint64_t)(bytes) << 0)
#define MATINT_Y_REG(reg) (MATINT_Y_OFFSET((reg)*64))

#define MATINT_Z_ROW(row) ((uint64_t)(row) << 20)

#define MATINT_ENABLE_MASK_IS_X (0ULL << 25)
#define MATINT_ENABLE_MASK_IS_Y (1ULL << 25)

#define MATINT_X_Y_ENABLE_VALUE(value) ((uint64_t)(value) << 32)

#define MATINT_X_Y_ENABLE_MODE(mode) ((uint64_t)(mode) << 38)

// genlut

#define GENLUT_SRC_OFFSET(bytes) ((uint64_t)(bytes) << 0)
#define GENLUT_SRC_REG(reg) (GENLUT_SRC_OFFSET((reg)*64))

#define GENLUT_SRC_X (0ULL << 10)
#define GENLUT_SRC_Y (1ULL << 10)

#define GENLUT_DST_REG(reg) ((uint64_t)(reg) << 20)

#define GENLUT_DST_X (0ULL << 25)
#define GENLUT_DST_Y (1ULL << 25)
#define GENLUT_DST_Z (1ULL << 26)

#define GENLUT_MODE_GENERATE_I16_U5 (4ULL << 53)
#define GENLUT_MODE_GENERATE_U16_U5 (6ULL << 53)
#define GENLUT_MODE_LOOKUP_16_U5 (14ULL << 53)

#define GENLUT_TABLE_X (0ULL << 59)
#define GENLUT_TABLE_Y (1ULL << 59)

#define GENLUT_TABLE_REG(reg) ((uint64_t)(reg) << 60)

#endif  // AMX_H
