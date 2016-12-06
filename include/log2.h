#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#include "useful_macros.h"

inline int get_log2_upper_power2(uint32_t pos_rounded_power_of_2);
inline uint32_t round_up_next_highest_power2(uint32_t pos);
inline int get_nb_bytes_power2_comp(uint32_t pos);
inline int get_nb_bytes_power2_annot(uint32_t pos);
inline int get_nb_bytes_power2_annot_bis(uint32_t pos, uint32_t pow2_pos);

inline int get_log2_upper_power2(uint32_t pos_rounded_power_of_2){

    const int MultiplyDeBruijnBitPosition_power2[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };

    return MultiplyDeBruijnBitPosition_power2[(uint32_t)(pos_rounded_power_of_2 * 0x077CB531U) >> 27];
}

inline uint32_t round_up_next_highest_power2(uint32_t pos){

    pos--;
    pos |= pos >> 1;
    pos |= pos >> 2;
    pos |= pos >> 4;
    pos |= pos >> 8;
    pos |= pos >> 16;
    pos++;

    return pos;
}

inline int get_nb_bytes_power2_comp(uint32_t pos){

    const int MultiplyDeBruijnBitPosition_power2[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };

    int log2_upper_pos = MultiplyDeBruijnBitPosition_power2[(uint32_t)(pos * 0x077CB531U) >> 27];

    if (log2_upper_pos > 6) return CEIL(log2_upper_pos-6, 7) + 1;
    else return 1;
}

inline int get_nb_bytes_power2_annot(uint32_t pos){

    const int MultiplyDeBruijnBitPosition_power2[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };

    uint32_t pow2_pos = round_up_next_highest_power2(pos);

    return CEIL(MultiplyDeBruijnBitPosition_power2[(uint32_t)(pow2_pos * 0x077CB531U) >> 27] + (pow2_pos == pos), 6);
}

inline int get_nb_bytes_power2_annot_bis(uint32_t pos, uint32_t pow2_pos){

    const int MultiplyDeBruijnBitPosition_power2[32] =
    {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };

    return CEIL(MultiplyDeBruijnBitPosition_power2[(uint32_t)(pow2_pos * 0x077CB531U) >> 27] + (pow2_pos == pos), 6);
}
