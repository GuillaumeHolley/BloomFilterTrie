#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "useful_macros.h"

extern const uint8_t POPCOUNT_8bit[256];
extern const uint8_t rev[256];
extern const uint8_t rev_MSB[16];
extern const uint8_t rev_LSB[16];

inline uint8_t reverse_word_8(uint8_t v);
inline int popcnt_8(uint8_t v);
inline int popcnt_8_par(const uint8_t* v, int start, int end);

inline uint8_t reverse_word_8(uint8_t v){

    return rev[v];
}

inline int popcnt_8(uint8_t v){

    return (int) POPCOUNT_8bit[v];
}

inline int popcnt_8_par(const uint8_t* v, int start, int end){

    ASSERT_NULL_PTR(v,"popcnt_8_par()")

    int hamming_weight = 0;
    for (int i=start; i<end; i++) hamming_weight += (int) POPCOUNT_8bit[v[i]];

    return hamming_weight;
}


