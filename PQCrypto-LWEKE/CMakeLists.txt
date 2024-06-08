foreach(PARAM_SET ${FRODOKEM_PARAMETER_SETS})
    set(SOURCES_FRODOKEM_${PARAM_SET}_SHAKE
        FrodoKEM/src/util.c FrodoKEM/src/frodo${PARAM_SET}.c FrodoKEM/src/frodo${PARAM_SET}_4x.c)
    set(SOURCES_FRODOKEM_${PARAM_SET}_AES ${SOURCES_FRODOKEM_${PARAM_SET}_SHAKE} ${CMAKE_SOURCE_DIR}/aes/aes_armv8.c)
endforeach()

set(SOURCES_FRODOKEM_ref "")
set(SOURCES_FRODOKEM_opt "")
set(SOURCES_FRODOKEM_neon "")
set(SOURCES_FRODOKEM_opt_amx
    ${CMAKE_SOURCE_DIR}/amx/frodokem_opt_matmul.c ${CMAKE_SOURCE_DIR}/amx/frodokem_opt_matmul_4x.c
    ${CMAKE_SOURCE_DIR}/amx/frodokem_matmul_4x.c)

set(FRODOKEM_VARIANTS AES SHAKE)
set(FRODOKEM_IMPLS ref opt neon opt_amx)
set(PQCTESTKAT_SUFFIXES "" _shake)
set(A_GENERATORS _AES128_FOR_A_ _SHAKE128_FOR_A_)

macro(ADD_KAT_TESTS_FRODOKEM)
    add_test(
        NAME ${LIBRARY}.KATs_match_spec
        COMMAND ${PQCGENKAT_KEM}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/PQCrypto-LWEKE/FrodoKEM)
endmacro()

foreach(PARAM_SET KAT_NUM IN ZIP_LISTS FRODOKEM_PARAMETER_SETS FRODOKEM_KAT_NUMS)
    foreach(VARIANT PQCTESTKAT_SUFFIX A_GENERATOR IN ZIP_LISTS FRODOKEM_VARIANTS PQCTESTKAT_SUFFIXES A_GENERATORS)
        set(PARAMETER_SET ${PARAM_SET}_${VARIANT})
        set(PQCTESTKAT_FILE PQCtestKAT_kem${PARAM_SET}${PQCTESTKAT_SUFFIX}.c)

        foreach(ALLOC stack mmap)
            foreach(IMPL ${FRODOKEM_IMPLS})
                set(LIBRARY frodokem_${PARAMETER_SET}_${ALLOC}_ref_${IMPL})

                set(SPEED_PREFIX speed)
                set(SPEED4X_PREFIX speed4x)

                # set(SPEED_POLYMUL_PREFIX speed_polymul)
                set(PQCGENKAT_KEM PQCgenKAT_kem_${LIBRARY})
                set(TEST test_${LIBRARY})

                add_library(${LIBRARY} STATIC ${HASH_SOURCES})
                target_compile_definitions(${LIBRARY} PUBLIC NIX _ARM_ _FAST_GENERIC_ ${A_GENERATOR})

                if(ALLOC STREQUAL stack)
                    target_compile_definitions(${LIBRARY} PUBLIC ALLOC_STACK)
                else()
                    target_compile_definitions(${LIBRARY} PUBLIC ALLOC_MMAP)
                endif()

                if(IMPL STREQUAL neon)
                    target_compile_definitions(${LIBRARY} PUBLIC OPT_NEON)
                endif()

                if(IMPL MATCHES amx)
                    target_compile_definitions(${LIBRARY} PUBLIC AMX)
                endif()

                target_include_directories(${LIBRARY} PUBLIC ${AMX_PATH} ${HASH_PATH})

                if(IMPL MATCHES ref)
                    target_compile_definitions(${LIBRARY} PUBLIC USE_REFERENCE)
                endif()

                target_link_libraries(${LIBRARY} PUBLIC neon_rng)

                foreach(SOURCE ${SOURCES_FRODOKEM_${PARAMETER_SET}} ${SOURCES_FRODOKEM_${IMPL}})
                    target_sources(${LIBRARY} PRIVATE ${SOURCE})
                endforeach()

                target_include_directories(${LIBRARY} PUBLIC FrodoKEM/src ${HASH_PATH} ${RAND_PATH})

                # foreach(SPEED_PREFIX ${SPEED_PREFIXES})
                set(SPEED ${SPEED_PREFIX}_${LIBRARY})
                add_executable_with_symlink(${SPEED} ${CMAKE_SOURCE_DIR}/speed/speed_kem.c)
                target_include_directories(${SPEED} PRIVATE ${RAND_PATH})
                target_compile_options(${SPEED} PRIVATE
                    -DNTESTS=100 --include=api_frodo${PARAM_SET}.h
                    -Dcrypto_kem_keypair=crypto_kem_keypair_Frodo${PARAM_SET}
                    -Dcrypto_kem_enc=crypto_kem_enc_Frodo${PARAM_SET}
                    -Dcrypto_kem_dec=crypto_kem_dec_Frodo${PARAM_SET})
                target_link_libraries(${SPEED} PRIVATE ${LIBRARY} cycles)

                set(SPEED4X ${SPEED4X_PREFIX}_${LIBRARY})
                add_executable_with_symlink(${SPEED4X} ${CMAKE_SOURCE_DIR}/speed/speed_kem4x.c)
                target_include_directories(${SPEED4X} PRIVATE ${RAND_PATH})
                target_compile_options(${SPEED4X} PRIVATE
                    -DNTESTS=100 --include=api_frodo${PARAM_SET}.h --include=api_frodo${PARAM_SET}_4x.h
                    -Dcrypto_kem_keypair=crypto_kem_keypair_Frodo${PARAM_SET}
                    -Dcrypto_kem_enc4x=crypto_kem_enc4x_Frodo${PARAM_SET}
                    -Dcrypto_kem_dec4x=crypto_kem_dec4x_Frodo${PARAM_SET})
                target_link_libraries(${SPEED4X} PRIVATE ${LIBRARY} cycles)

                add_executable_with_symlink(${TEST}
                    ${CMAKE_SOURCE_DIR}/test/test_kem.cpp ${CMAKE_SOURCE_DIR}/test/test_kem4x.cpp)
                target_include_directories(${TEST} PRIVATE ${RAND_PATH})
                target_compile_options(${TEST} PRIVATE
                    -DTEST_NAME=${LIBRARY} --include=api_frodo${PARAM_SET}.h --include=api_frodo${PARAM_SET}_4x.h
                    -Dcrypto_kem_keypair=crypto_kem_keypair_Frodo${PARAM_SET}
                    -Dcrypto_kem_enc=crypto_kem_enc_Frodo${PARAM_SET}
                    -Dcrypto_kem_dec=crypto_kem_dec_Frodo${PARAM_SET}
                    -Dcrypto_kem_enc4x=crypto_kem_enc4x_Frodo${PARAM_SET}
                    -Dcrypto_kem_dec4x=crypto_kem_dec4x_Frodo${PARAM_SET})
                target_link_libraries(${TEST} PRIVATE ${LIBRARY} GTest::gtest_main)
                gtest_discover_tests(${TEST} DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})

                add_executable_with_symlink(${PQCGENKAT_KEM} FrodoKEM/tests/${PQCTESTKAT_FILE})
                target_compile_options(${PQCGENKAT_KEM} PRIVATE -Wno-unused-result)
                target_link_libraries(${PQCGENKAT_KEM} PRIVATE ${LIBRARY})

                ADD_KAT_TESTS_FRODOKEM()
            endforeach()
        endforeach()
    endforeach()
endforeach()

set(SAMPLE_IMPLS opt neon opt_amx)

foreach(PARAMS_N 640 976 1344)
    add_executable_with_symlink(test_frodokem_${PARAMS_N}_sample ${CMAKE_SOURCE_DIR}/test/test_sample.cpp)
    target_compile_definitions(test_frodokem_${PARAMS_N}_sample PRIVATE PARAMS_N=${PARAMS_N} ALLOC_MMAP)
    target_include_directories(test_frodokem_${PARAMS_N}_sample PRIVATE ${AMX_PATH} test FrodoKEM/src)
    target_link_libraries(test_frodokem_${PARAMS_N}_sample PRIVATE GTest::gtest_main)
    gtest_discover_tests(test_frodokem_${PARAMS_N}_sample DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})

    set(SPEED_SAMPLE speed_sample_frodokem_${PARAMS_N})
    add_executable_with_symlink(${SPEED_SAMPLE} ${CMAKE_SOURCE_DIR}/speed/speed_sample.c)
    target_include_directories(${SPEED_SAMPLE} PRIVATE ${AMX_PATH} ${RAND_PATH} FrodoKEM/src)
    target_compile_definitions(${SPEED_SAMPLE} PRIVATE ALLOC_MMAP NTESTS=1024 PARAMS_N=${PARAMS_N})
    target_link_libraries(${SPEED_SAMPLE} PRIVATE neon_rng cycles)

    foreach(SAMPLE_IMPL ${SAMPLE_IMPLS})
        set(SAMPLE_CT_EXPERIMENT sample_ct_experiment_frodokem_${PARAMS_N}_${SAMPLE_IMPL})
        add_executable_with_symlink(${SAMPLE_CT_EXPERIMENT} ${CMAKE_SOURCE_DIR}/speed/sample_ct_experiment.c)
        target_include_directories(${SAMPLE_CT_EXPERIMENT} PRIVATE ${AMX_PATH} ${RAND_PATH} FrodoKEM/src)
        target_compile_definitions(${SAMPLE_CT_EXPERIMENT} PRIVATE ALLOC_MMAP NTESTS=1024 PARAMS_N=${PARAMS_N})

        if(SAMPLE_IMPL STREQUAL neon)
            target_compile_definitions(${SAMPLE_CT_EXPERIMENT} PRIVATE OPT_NEON)
        elseif(SAMPLE_IMPL STREQUAL opt_amx)
            target_compile_definitions(${SAMPLE_CT_EXPERIMENT} PRIVATE AMX)
        endif()

        target_link_libraries(${SAMPLE_CT_EXPERIMENT} PRIVATE neon_rng cycles)
    endforeach()

    foreach(VARIANT A_GENERATOR IN ZIP_LISTS FRODOKEM_VARIANTS A_GENERATORS)
        foreach(ALLOC stack mmap)
            foreach(IMPL ${FRODOKEM_IMPLS})
                foreach(BATCHING "" 4x)
                    if(BATCHING STREQUAL 4x)
                        set(_BATCHING _4x)
                    else()
                        set(_BATCHING "")
                    endif()

                    set(SPEED_MATMUL speed${BATCHING}_matmul_frodokem_${PARAMS_N}_${VARIANT}_${ALLOC}_ref_${IMPL})

                    if(IMPL STREQUAL opt OR IMPL STREQUAL neon)
                        add_executable_with_symlink(${SPEED_MATMUL}
                            ${CMAKE_SOURCE_DIR}/speed/speed_frodokem_opt_neon_matmul${_BATCHING}.c)
                    else()
                        add_executable_with_symlink(${SPEED_MATMUL}
                            ${CMAKE_SOURCE_DIR}/speed/speed_frodokem_${IMPL}_matmul${_BATCHING}.c)
                    endif()

                    target_include_directories(${SPEED_MATMUL} PRIVATE ${AMX_PATH} ${RAND_PATH} FrodoKEM/src)
                    target_compile_definitions(${SPEED_MATMUL} PRIVATE ${A_GENERATOR} NTESTS=1024 PARAMS_N=${PARAMS_N})

                    if(IMPL STREQUAL neon)
                        target_compile_definitions(${SPEED_MATMUL} PRIVATE OPT_NEON)
                    endif()

                    if(ALLOC STREQUAL mmap)
                        target_compile_definitions(${SPEED_MATMUL} PRIVATE ALLOC_MMAP)
                    else()
                        target_compile_definitions(${SPEED_MATMUL} PRIVATE ALLOC_STACK)
                    endif()

                    target_link_libraries(${SPEED_MATMUL} PRIVATE
                        neon_rng cycles frodokem_${PARAMS_N}_${VARIANT}_${ALLOC}_ref_${IMPL})
                endforeach()
            endforeach()
        endforeach()
    endforeach()
endforeach()
