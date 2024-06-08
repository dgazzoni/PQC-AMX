cmake_minimum_required(VERSION 3.18)

project(PQC-AMX C CXX ASM)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

option(CMAKE_UNITY_BUILD "merge source files into batches for faster compilation" ON)
option(SANITIZER "build with AddressSanitizer and UBSanitizer support" OFF)
option(BUILD_TESTING "build with tests enabled" ON)

set(CMAKE_UNITY_BUILD_BATCH_SIZE 0)

if(BUILD_TESTING)
    include(CTest)
endif()

include(BuildCaching.cmake)
include(CreateSymlink.cmake)
include(CPUDetection.cmake)

add_compile_options(-Wall -Wextra -pedantic)

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-language-extension-token -Wno-strict-prototypes)
endif()

if((CMAKE_C_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 15) OR
    CMAKE_C_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wno-array-parameter)
endif()

if(SANITIZER)
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

add_compile_options(-march=armv8-a+crypto)

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_compile_options(
        "SHELL:-mllvm -align-all-functions=6" "SHELL:-mllvm -align-all-nofallthru-blocks=6")
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")
    add_compile_options(-falign-functions=64)
endif()

if(APPLE)
    execute_process(COMMAND brew --prefix openssl OUTPUT_VARIABLE OPENSSL_ROOT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

find_package(OpenSSL REQUIRED)

set(BUILD_GMOCK OFF CACHE BOOL "")
set(INSTALL_GTEST OFF CACHE BOOL "")
add_subdirectory(googletest)

include(GoogleTest)
set(GTEST_DISCOVERY_TIMEOUT 60)

add_library(ref_rng OBJECT neon-ntt/saber/m1_benchmarks/rng.c)
target_compile_definitions(ref_rng PUBLIC randombytes_init=nist_randombytes_init randombytes=nist_randombytes)
target_link_libraries(ref_rng PRIVATE OpenSSL::SSL OpenSSL::Crypto)
target_compile_options(ref_rng PRIVATE -Wno-sign-compare -Wno-unused-parameter)

if((CMAKE_C_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER 10) OR
    CMAKE_C_COMPILER_ID MATCHES "GNU")
    target_compile_options(ref_rng PRIVATE -Wno-unused-but-set-variable)
endif()

# target_link_libraries(ref_rng PUBLIC OpenSSL::SSL OpenSSL::Crypto)
if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_library(opt_rng OBJECT rng_opt/rng_inline_asm.c)
else()
    add_library(opt_rng OBJECT rng_opt/rng.c)
    set_source_files_properties(rng_opt/rng.c PROPERTIES COMPILE_FLAGS -fno-strict-aliasing)
endif()

target_compile_definitions(opt_rng PUBLIC randombytes_init=opt_randombytes_init randombytes=opt_randombytes)
add_library(neon_rng ALIAS opt_rng)

set(FRODOKEM_PARAMETER_SETS 640 976 1344)
set(FRODOKEM_KAT_NUMS 19888 31296 43088)

set(SABER_PARAMETER_SETS lightsaber saber firesaber)
set(SABER_KAT_NUMS 1568 2304 3040)

set(HASH_PATH ${CMAKE_SOURCE_DIR}/neon-ntt/common)
set(HASH_SOURCES ${HASH_PATH}/fips202.c ${HASH_PATH}/fips202x2.c ${HASH_PATH}/feat.S)

set(RAND_PATH ${CMAKE_SOURCE_DIR}/rng_opt)

set(AMX_PATH ${CMAKE_SOURCE_DIR}/amx)

add_library(cycles OBJECT ${CMAKE_SOURCE_DIR}/neon-ntt/common/m1cycles.c)

target_include_directories(cycles PUBLIC ${CMAKE_SOURCE_DIR}/neon-ntt/common)

add_subdirectory(neon-ntt)
add_subdirectory(PQCrypto-LWEKE)

# Tests
add_executable(test_rng test/test_rng.cpp)
target_include_directories(test_rng PRIVATE reference/Reference_Implementation/crypto_kem/ntruhps2048509)

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    target_compile_options(test_rng PRIVATE -Wno-macro-redefined)
endif()

target_link_libraries(test_rng PRIVATE ref_rng opt_rng GTest::gtest_main)
gtest_discover_tests(test_rng DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})

add_library(aes_armv8_for_test OBJECT aes/aes_armv8.c)
target_compile_definitions(aes_armv8_for_test PRIVATE
    AES128_load_schedule=ARMv8_AES128_load_schedule AES128_ECB_enc_sch=ARMv8_AES128_ECB_enc_sch
    AES128_free_schedule=ARMv8_AES128_free_schedule)

add_executable(test_aes test/test_aes.cpp PQCrypto-LWEKE/common/aes/aes_c.c)
target_include_directories(test_aes PRIVATE PQCrypto-LWEKE/common/aes aes)
target_link_libraries(test_aes PRIVATE aes_armv8_for_test GTest::gtest_main)
gtest_discover_tests(test_aes DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})

foreach(TEST_TYPE polymodmul polymodmul_matmul)
    add_executable(test_${TEST_TYPE} test/test_polymodmul.cpp amx/aux_routines.c amx/${TEST_TYPE}.c)
    target_include_directories(test_${TEST_TYPE} PRIVATE ${AMX_PATH} test)
    target_link_libraries(test_${TEST_TYPE} PRIVATE GTest::gtest_main)
    gtest_discover_tests(test_${TEST_TYPE} DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})
endforeach()

target_compile_definitions(test_polymodmul_matmul PRIVATE POLYMODMUL_MATMUL)

# Benchmarks
add_executable(speed_aes speed/speed_aes.c aes/aes_armv8.c)
target_include_directories(speed_aes PRIVATE PQCrypto-LWEKE/common/aes)
target_compile_definitions(speed_aes PRIVATE NTESTS=1024)
target_link_libraries(speed_aes PRIVATE cycles)
