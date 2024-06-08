set(SABER_SOURCES_COMMON SABER_indcpa.c kem.c pack_unpack.c verify.c cbd.c __asm_pack_unpack.S)
set(SABER_SOURCES_neon __asm_NTT.S __asm_iNTT.S __asm_mul.S __asm_narrow.S)
set(SABER_SOURCES_amx_polymul poly.c poly_matrixvectormulround.c)
set(SABER_SOURCES_amx_matmul poly.c)

set(AMX_SOURCES_amx_polymul ${CMAKE_SOURCE_DIR}/amx/aux_routines.c ${CMAKE_SOURCE_DIR}/amx/polymodmul.c)

# We use the polymul code for single polynomial multiplication in the matmul implementation; matmul only reimplements
# the matrix-vector product routine with a faster algorithm.
set(AMX_SOURCES_amx_matmul ${CMAKE_SOURCE_DIR}/amx/aux_routines.c ${CMAKE_SOURCE_DIR}/amx/polymodmul.c)
set(AMX_SOURCES_amx_matmul_lightsaber ${CMAKE_SOURCE_DIR}/amx/poly_lightsaber_amx.c)
set(AMX_SOURCES_amx_matmul_saber ${CMAKE_SOURCE_DIR}/amx/poly_saber_amx.c)
set(AMX_SOURCES_amx_matmul_firesaber ${CMAKE_SOURCE_DIR}/amx/poly_firesaber_amx.c)

macro(ADD_KAT_TESTS_SABER)
    # https://stackoverflow.com/a/3071370/523079
    add_test(
        NAME ${LIBRARY}.KATs_match_spec
        COMMAND
        ${CMAKE_COMMAND}
        -DKATgen_cmd=${PQCGENKAT_KEM}
        -DKAT_expected=${CMAKE_SOURCE_DIR}/neon-ntt/${PARAMETER_SET}/KAT/PQCkemKAT_${KAT_NUM}
        -DKAT_actual=PQCkemKAT_${KAT_NUM}
        -DWORKING_DIRECTORY=${CMAKE_BINARY_DIR}/KAT/${LIBRARY}
        -DBUILD_DIRECTORY=${CMAKE_BINARY_DIR}
        -DSRC_DIRECTORY=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_SOURCE_DIR}/CompareKATs.cmake)
endmacro()

set(IMPLS neon amx_polymul amx_matmul)
set(IMPL_DIRS m1_benchmarks amx amx)

foreach(PARAMETER_SET KAT_NUM IN ZIP_LISTS SABER_PARAMETER_SETS SABER_KAT_NUMS)
    foreach(ALLOC stack mmap)
        foreach(IMPL IMPL_DIR IN ZIP_LISTS IMPLS IMPL_DIRS)
            set(LIBRARY ${PARAMETER_SET}_${ALLOC}_BHK_${IMPL})

            set(SPEED_PREFIX speed)
            set(SPEED_MATRIXVECTORMULROUND_PREFIX speed_matrixvectormulround)

            set(PQCGENKAT_KEM PQCgenKAT_kem_${LIBRARY})
            set(TEST test_${LIBRARY})

            add_library(${LIBRARY} STATIC ${HASH_SOURCES})

            target_include_directories(${LIBRARY} PUBLIC ${PARAMETER_SET}/${IMPL_DIR})
            target_link_libraries(${LIBRARY} PUBLIC neon_rng)

            foreach(SOURCE ${SABER_SOURCES_COMMON} ${SABER_SOURCES_${IMPL}})
                target_sources(${LIBRARY} PRIVATE ${PARAMETER_SET}/${IMPL_DIR}/${SOURCE})
            endforeach()

            if(ALLOC STREQUAL stack)
                target_compile_definitions(${LIBRARY} PUBLIC ALLOC_STACK)
            else()
                target_compile_definitions(${LIBRARY} PUBLIC ALLOC_MMAP)
            endif()

            if(IMPL MATCHES amx)
                target_sources(${LIBRARY} PRIVATE ${AMX_SOURCES_${IMPL}} ${AMX_SOURCES_${IMPL}_${PARAMETER_SET}})
            endif()

            target_include_directories(${LIBRARY} PUBLIC
                ${PARAMETER_SET}/${IMPL_DIR} ${AMX_PATH} ${HASH_PATH} ${RAND_PATH})

            set(SPEED ${SPEED_PREFIX}_${LIBRARY})
            add_executable_with_symlink(${SPEED} ${CMAKE_SOURCE_DIR}/speed/speed_kem.c)
            target_compile_options(${SPEED} PRIVATE -DNTESTS=1024 --include=api.h)
            target_include_directories(${SPEED} PRIVATE ${RAND_PATH})
            target_link_libraries(${SPEED} PRIVATE ${LIBRARY} cycles)

            set(SPEED_MATRIXVECTORMULROUND ${SPEED_MATRIXVECTORMULROUND_PREFIX}_${LIBRARY})

            if(IMPL MATCHES neon)
                add_executable_with_symlink(${SPEED_MATRIXVECTORMULROUND}
                    ${CMAKE_SOURCE_DIR}/speed/speed_saber_matrixvectormulround_neon.c)
            else()
                add_executable_with_symlink(${SPEED_MATRIXVECTORMULROUND}
                    ${CMAKE_SOURCE_DIR}/speed/speed_saber_matrixvectormulround_amx.c)
            endif()

            target_compile_options(${SPEED_MATRIXVECTORMULROUND} PRIVATE -DNTESTS=1024)
            target_include_directories(${SPEED_MATRIXVECTORMULROUND} PRIVATE ${RAND_PATH})
            target_link_libraries(${SPEED_MATRIXVECTORMULROUND} PRIVATE ${LIBRARY} cycles)

            add_executable_with_symlink(${TEST} ${CMAKE_SOURCE_DIR}/test/test_kem.cpp)
            target_include_directories(${TEST} PRIVATE ${RAND_PATH})
            target_compile_options(${TEST} PRIVATE -DTEST_NAME=${LIBRARY} --include=api.h)
            target_link_libraries(${TEST} PRIVATE ${LIBRARY} GTest::gtest_main)
            gtest_discover_tests(${TEST} DISCOVERY_TIMEOUT ${GTEST_DISCOVERY_TIMEOUT})

            add_executable_with_symlink(${PQCGENKAT_KEM} common/PQCgenKAT_kem.c)
            target_compile_options(${PQCGENKAT_KEM} PRIVATE -Wno-unused-result)
            target_link_libraries(${PQCGENKAT_KEM} PRIVATE ${LIBRARY})

            ADD_KAT_TESTS_SABER()
        endforeach()
    endforeach()
endforeach()

set(SPEED_POLYMODMUL speed_polymodmul_saber_amx)
add_executable_with_symlink(${SPEED_POLYMODMUL}
    ${CMAKE_SOURCE_DIR}/speed/speed_saber_polymodmul_amx.c ${CMAKE_SOURCE_DIR}/amx/aux_routines.c
    ${CMAKE_SOURCE_DIR}/amx/polymodmul.c)
target_include_directories(${SPEED_POLYMODMUL} PRIVATE ${AMX_PATH} ${RAND_PATH})
target_compile_definitions(${SPEED_POLYMODMUL} PRIVATE NTESTS=1024)
target_link_libraries(${SPEED_POLYMODMUL} PRIVATE neon_rng cycles)

set(SPEED_POLYMODMUL_MATMUL speed_polymodmul_matmul_saber_amx)
add_executable_with_symlink(${SPEED_POLYMODMUL_MATMUL}
    ${CMAKE_SOURCE_DIR}/speed/speed_saber_polymodmul_amx.c ${CMAKE_SOURCE_DIR}/amx/polymodmul_matmul.c)
target_include_directories(${SPEED_POLYMODMUL_MATMUL} PRIVATE ${AMX_PATH} ${RAND_PATH})
target_compile_definitions(${SPEED_POLYMODMUL_MATMUL} PRIVATE ALLOC_MMAP NTESTS=1024)
target_link_libraries(${SPEED_POLYMODMUL_MATMUL} PRIVATE neon_rng cycles)
