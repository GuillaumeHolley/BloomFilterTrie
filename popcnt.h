#ifndef DEF_POPCNT
#define DEF_POPCNT

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "useful_macros.h"

inline uint8_t reverse_word_8(uint8_t v) {
    const uint8_t rev_MSB[16] = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
    const uint8_t rev_LSB[16] = {0, 64, 128, 192, 16, 80, 144, 208, 32, 96, 160, 224, 48, 112, 176, 240};
    return rev_LSB[v & 0xf] | rev_MSB[v >> 4];
}

inline int popcnt_8(uint8_t v) {
    const uint8_t tbl[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
    return tbl[v & 0xf] + tbl[v >> 4];
}

inline int popcnt_8_par(uint8_t* v, int start, int end) {

    ASSERT_NULL_PTR(v,"popcnt_8_par()")

    int hamming_weight = 0;
    int i=0;
    const uint8_t tbl[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

    for (i=start; i<end; i++) hamming_weight += tbl[v[i] & 0xf] + tbl[v[i] >> 4];

    return hamming_weight;
}

#endif


