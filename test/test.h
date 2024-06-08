#ifndef TEST_H
#define TEST_H

#include <cstddef>
#include <cstdlib>

#include "gtest/gtest.h"

// https://stackoverflow.com/a/10062016/523079
template <typename T, size_t size>
::testing::AssertionResult ArraysMatch(const T (&expected)[size], const T (&actual)[size]) {
    for (size_t i(0); i < size; ++i) {
        if (expected[i] != actual[i]) {
            return ::testing::AssertionFailure()
                   << "expected[" << i << "] (" << expected[i] << ") != actual[" << i << "] (" << actual[i] << ")";
        }
    }

    return ::testing::AssertionSuccess();
}

template <class T1, class T2>
::testing::AssertionResult ArraysMatch(const T1 *expected, const T2 *actual, size_t array_len) {
    for (size_t i(0); i < array_len; ++i) {
        if (expected[i] != actual[i]) {
            ::testing::Message msg;

            msg << "expected[" << i << "] (" << +expected[i] << ") != actual[" << i << "] (" << +actual[i] << ")";
            return ::testing::AssertionFailure() << msg;
        }
    }

    return ::testing::AssertionSuccess();
}

template <class T1, class T2>
::testing::AssertionResult ArraysMatchAllDifferences(const T1 *expected, const T2 *actual, size_t array_len) {
    bool error = false;
    ::testing::Message msg;
    msg << "Errors: ";
    for (size_t i(0); i < array_len; ++i) {
        if (expected[i] != actual[i]) {
            msg << i << ", ";

            error = true;
        }
    }

    if (error) {
        return ::testing::AssertionFailure() << msg;
    }
    else {
        return ::testing::AssertionSuccess();
    }
}

static inline void gen_random(uint16_t v[], size_t length) {
    for (size_t i = 0; i < length; i++) {
        v[i] = rand() & 0xFFFF;
    }
}

#endif  // TEST_H
